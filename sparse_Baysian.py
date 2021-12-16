# -*- coding: utf-8 -*-
"""
Created on Mon Dec 13 21:16:00 2021

@author: huyanfeng
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import random as rd
import copy
import time
import math
import heapq

def get_distance(K,radius):
    '''
    this function is used to get the distances between users and BS
    Parameters
    ----------
    K : Int
        the number of total users.
    radius : double
        the radius of a cell.

    Returns
    -------
    distance : nparray of K*1,double
        the distance between users and BS.
    '''
    if(K<=0):
        return np.zeros((0,0))
    else:
        distance = radius*np.abs(np.random.rand(K,1))
        return distance
    

def get_Channel(p,distance,d0,pl_exp,K,M):
    '''
    Parameters
    ----------
    p : double, unit of dbm
        power of transmit signals.
    distance : nparray
        distance of users and BS.
    d0 : double
        reference distance.
    pl_exp : double
        path loss factor.
    K : int
        number of total users.
    M : int
        number of antennas.

    Returns
    -------
    G : nparray
        a real channel matrix.
    '''
    c0 = 34.5 #+ 20*np.log10(d0)      #总感觉这个参考点的取值是不是有点大
    noise_figure = 9
    n0 = -174
    bw = 10e6
    snr = p-c0-noise_figure-10*np.log10(bw)-n0
    snr = 10**(snr/10.0)
    pow_lsf = 2*snr*(1/(1+distance/d0)**pl_exp)
    h = (np.random.randn(M,K)+1j*np.random.randn(M,K))/np.sqrt(2)
    G = np.transpose(np.sqrt(pow_lsf).repeat(M,axis=1))*h
    return np.transpose(G)

def get_active_user(K,Ka):
    ''''
    获取活跃用户索引
    '''
    #rd.seed(0)
    user_index = rd.sample(range(K), Ka)
    user_index.sort()
    return user_index

def get_codeBook(L,N):
    '''
    获取独立同分布的公共码本
    L为压缩后的向量尺寸，N为压缩前的向量尺寸
    '''
    #np.random.seed(0)
    A = (np.random.randn(L,N)+1j*np.random.randn(L,N))/np.sqrt(2)
    return A

def get_sparseVector(user_index,Ka,K,N):
    '''活跃用户产生稀疏向量'''
    v = np.zeros((N,K))
    for i in range(Ka):
        v[rd.randint(0, N-1),user_index[i]]=1.0
    return v

def sparse_baysian_learning(A,y,iterNum,alpha0,a,b):
    '''
    该函数用以利用稀疏贝叶斯学习方法估计活跃用户信道矩阵
    Parameters
    ----------
    A : nparray
        这是感知矩阵，收发端均已知.
    y : nparray
        接收信号矩阵.
    iterNum : int
        迭代次数.
    alpha0 : double
        噪声方差的倒数.
    a : double
        gamma分布第一参数.
    b : double
        gamma分布第二参数.

    Returns
    -------
    u : nparry
        估计信道均值矩阵.
    alpha : nparry
        估计信道方差倒数向量.
    sigma : nparry
        协方差矩阵.
    '''
    (L,N) = np.shape(A)
    M = np.shape(y)[1]
    u = np.zeros((N,M))
    alpha = np.ones((N,1))
    for i in range(iterNum):
        B = np.eye(N)*alpha
        sigma = np.linalg.inv(alpha0*np.dot(A.conj().T,A)+B)
        u = np.dot(np.dot(alpha0*sigma,A.conj().T),y)
        alpha = (1+2*a)/((np.diag(sigma)).reshape(N,1)+((u*u.conj()).sum(axis=1)/M).reshape(N,1)+2*b)
    return u,alpha,sigma

def channel_estimation(A,user_index,K,Ka,M,N,G,a,b,iterNum,repeatNum,alpha0):
    '''
    信道估计函数，在原先稀疏贝叶斯学习的基础上，重复随机发送压缩感知向量，
    从中选取出现次数最多的前Ka个信道向量作为活跃用户信道向量。
    这里的出现次数指的是误差不超过5%的估计向量，表示估计结果一致
    输入参数含义与sparse_baysian_learning相一致
    输出参数为估计得到的信道矩阵
    '''
    u = np.zeros((repeatNum,Ka,M))+1j*np.zeros((repeatNum,Ka,M))
    channel_list = []
    estimation_Num = []
    for i in range(repeatNum):
        x = get_sparseVector(user_index, Ka, K, N)    #每次重复发送都随机生成稀疏向量
        #加上方差为1的噪声向量
        y = np.dot(np.dot(A,x),G)+(np.random.randn(L,M)+1j*np.random.randn(L,M))/np.sqrt(2) 
        ux,alpha,sigma = sparse_baysian_learning(A, y, iterNum, alpha0, a, b)
        #选出方差前Ka小的估计结果作为稀疏贝叶斯学习的估计信道值
        temp = heapq.nsmallest(Ka,range(N),alpha.take)
        u[i,:,:] = ux[temp,:]
        for j in range(Ka):
            flag = 0
            #遍历信道列表，若误差在5%以内则放入同一信道列表之中
            for k in range(len(channel_list)):
                if(sum(np.abs(u[i,j,:]-np.sum(channel_list[k],axis=0)/len(channel_list[k]))/np.abs(u[i,j,:]))/M<=0.1):
                    channel_list[k].append(u[i,j,:])
                    flag = 1
                    break
            #否则新增信道列表
            if(flag == 0):
                channel_list.append([u[i,j,:]])
    #以下for循环用以统计估计得到的相同信道出现的次数，并求相同信道向量的均值作为信道向量
    for i in range(len(channel_list)):
        estimation_Num.append(len(channel_list[i]))
        channel_list[i] = np.sum(channel_list[i],axis=0)/len(channel_list[i])
    Inf = 0
    res = []
    #以下for循环用以取出出现次数最多的前Ka个信道向量作为输出
    for i in range(Ka):
        res.append(channel_list[estimation_Num.index(max(estimation_Num))])
        estimation_Num[estimation_Num.index(max(estimation_Num))] = Inf
    return np.array(res)

def MSE(H,G):
    '''
    均方误差计算函数，输入参数为真实信道与估计信道，计算二者的均方误差与相对误差
    '''
    if(len(H)==0 or len(G)==0 or np.shape(H)!=np.shape(G)):
        return np.zeros((0,0))
    else:
        temp1 = (np.dot(H,H.conj().T)).diagonal()
        temp2 = (np.dot(G,G.conj().T)).diagonal()
        H = H[temp1.argsort()]
        G = G[temp2.argsort()]
        mse = np.average((H-G)*((H-G).conj()),axis=1)
        error_ratio = np.sqrt(mse)/(np.average(np.abs(G),axis=1))
        return mse, error_ratio

'''def VAMP(A,y,iterNum,mu,mu0):
    [M,Ka] = np.shape(A)
    L = np.shape(y)[1]
    '''
    


if __name__ == '__main__':
    '''
    K: 用户总数
    M: 天线数
    radius: 小区半径
    Ka: 活跃用户数
    L: 压缩后的向量长度
    N: 压缩前的向量长度
    p: 发送信号功率
    d0: 参考距离
    pl_exp: 路径损耗因子
    a: gamma分布参数1
    b: gamma分布参数2
    iterNum: 迭代循环次数
    alpha0: 噪声方差的倒数
    repeatNum: 重复发送次数
    '''
    K = 1024
    M = 128
    radius = 1000
    Ka = 24
    L = 100
    N = 512
    p = 17
    d0 = 10
    pl_exp = 3.7
    a = 1e-8
    b = 1e-8
    iterNum = 300
    alpha0 = 1
    repeatNum = 5
    #获取用户到基站距离
    distance = get_distance(K, radius)
    #获取真实信道向量
    G = get_Channel(p, distance, d0, pl_exp, K, M)
    #获取活跃用户下标索引
    user_index = get_active_user(K,Ka)
    #获取公共压缩感知码本
    A = get_codeBook(L, N)
    #进行信道估计过程
    print("start time:")
    print(time.strftime("%Y-%m-%d %H:%M:%S", time.localtime()))
    channels = channel_estimation(A, user_index, K, Ka, M, N, G, a, b, iterNum, repeatNum, alpha0)
    mse,error_ratio = MSE(channels,G[user_index,:])
    print("end time:")
    print(time.strftime("%Y-%m-%d %H:%M:%S", time.localtime()))
    # x = get_sparseVector(user_index, Ka, K, N)
    # y = np.dot(np.dot(A,x),G)+(np.random.randn(L,M)+1j*np.random.randn(L,M))/np.sqrt(2)
    # u,alpha,sigma = sparse_baysian_learning(A,y,iterNum,alpha0,a,b)
    # temp = G[user_index,:]
    # temp1 = heapq.nsmallest(Ka,range(N),alpha.take)
    # temp2 = u[temp1,:]
    plt.plot(np.array(range(Ka)),np.real(mse),label='MSE')
    plt.plot(np.array(range(Ka)),np.real(error_ratio),label='Error_ratio')
    plt.xlabel('user_index')
    plt.ylim([0,1])
    plt.legend()
    plt.show()
    
    
    
    
    