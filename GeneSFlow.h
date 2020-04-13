/*
 * @Author: Fater7 
 * @Date: 2019-07-04 16:16:58 
 * @Last Modified by: Fater7
 * @Last Modified time: 2019-07-04 23:15:41
 */
/*
 * 基于遗传算法的正函数极值计算框架
 */

#pragma once

#include <string>
#include <vector>
#include <cstdlib>
#include <ctime>
#include <cmath>
#include <algorithm>
#include <iostream>

#include "GeneBody.h"
#include "GSFRange.h"

using namespace std;

class GeneSFlow
{
public:
    GeneSFlow(int tM, int tT, GSFRange trange, double tPre, double tFx(double));
    void SetProbForCross(double d);
    void SetProbForMtp(double d);
    double start();

private:
    void CalIdvFitness(int groupid);
    void ChooseIdvToNext(int curGroupId, int nextGroupId);
    void CrossOver(int groupId);
    void IdvVariation(int groupid);
    double GSFDeCode(string genelist);
    string GetGeneWithLen(int length);
    int GSFRandom(int range);
    void PrintGroup(int groupid);
    void GSFLog(string log);


private:
    vector<GeneBody> GeneGroup[2];  //双种群，交替计算
    double (*Fx)(double x);         //适度函数
    GSFRange groupRange;            //种群范围
    int M;          //种群规模
    int T;          //进化代数
    double Pc;      //交叉概率
    double Pm;      //变异概率
    double Pre;     //目标精度
    int genelen;    //基因长度
    double genecount;  //genelen长度下的个体种类
};
