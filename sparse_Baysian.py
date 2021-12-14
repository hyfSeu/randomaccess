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
    rd.seed(0)
    user_index = rd.sample(range(K), Ka)
    user_index.sort()
    return user_index

def get_codeBook(L,N):
    '''
    获取独立同分布的公共码本
    L为压缩后的向量尺寸，N为压缩前的向量尺寸
    '''
    np.random.seed(0)
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
    
    


if __name__ == '__main__':
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
    distance = get_distance(K, radius)
    G = get_Channel(p, distance, d0, pl_exp, K, M)
    user_index = get_active_user(K,Ka)
    A = get_codeBook(L, N)
    x = get_sparseVector(user_index, Ka, K, N)
    y = np.dot(np.dot(A,x),G)+(np.random.randn(L,M)+1j*np.random.randn(L,M))/np.sqrt(2)
    u,alpha,sigma = sparse_baysian_learning(A,y,iterNum,alpha0,a,b)
    temp = G[user_index,:]
    temp1 = heapq.nsmallest(Ka,range(N),alpha.take)
    temp2 = u[temp1,:]
    
    
    
    
    