/*
 * @Author: Fater7 
 * @Date: 2019-07-04 16:17:04 
 * @Last Modified by: Fater7
 * @Last Modified time: 2019-07-04 23:16:16
 */

#include "GeneSFlow.h"

GeneSFlow::GeneSFlow(int tM, int tT, GSFRange trange, 
                     double tPre, double tFx(double))
{
    //初始化种群规模，进化代数，种群范围，适应度函数，目标精度
    M = 2 * tM;
    T = tT;
    groupRange = trange;
    Fx = tFx;
    Pre = tPre;

    //由精度与种群范围计算个体长度与个体种类
    double count = groupRange.length / Pre;
    genelen = 0;
    while(pow(2.0, genelen) < count)
    {
        genelen++;
    }
    genecount = pow(2, genelen);

    //初始化基因种群
    GeneGroup[0].clear();
    
    for(int i = 0; i < M; i++)
    {
        string tmpgene = GetGeneWithLen(genelen);
        GeneBody tmpGB(tmpgene);
        GeneGroup[0].push_back(tmpGB);
    }

    GSFLog("个体长度: " + to_string(genelen));
    GSFLog("初始种群:");
    PrintGroup(0);
}

/*
 * 运行遗传算法，返回最终结果
 */
double GeneSFlow::start()
{
    srand(time(NULL));

    for(int i = 0; i < T; i++)
    {
        GSFLog("Generation " + to_string(i+1) + ":");
        CalIdvFitness(i % 2);
        ChooseIdvToNext(i % 2, (i+1) % 2);
        GSFLog("适应度选择:");
        PrintGroup((i+1) % 2);

        CrossOver((i+1) % 2);
        GSFLog("交叉运算:");
        PrintGroup((i+1) % 2);

        IdvVariation((i+1) % 2);
        GSFLog("变异运算:");
        PrintGroup((i+1) % 2);
        GSFLog("\n");
    }

    CalIdvFitness(T % 2);
    return GSFDeCode(GeneGroup[T % 2][0].geneList);
}

/*
 * 按照交叉概率Pc对种群中的个体进行配对交叉
 */
void GeneSFlow::CrossOver(int gid)
{
    vector<int> idvUsed;
    for(int i = 0; i < M; i++)
    {
        if(count(idvUsed.begin(), idvUsed.end(), i) > 0)
        {
            //如果该个体已被配对，则跳过
            continue;
        }
        else
        {
            double tmpnum = GSFRandom(1000)/1000.0;
            if(tmpnum <= Pc)
            {
                //被概率选中，开始配对
                idvUsed.push_back(i);
                int gb2num;
                do
                {
                    gb2num = GSFRandom(M);
                } while (count(idvUsed.begin(), idvUsed.end(), gb2num) > 0);
                
                idvUsed.push_back(gb2num);

                int crosslen = GSFRandom(genelen-1)+1;
                string idv1part = GeneGroup[gid][i].geneList.substr(0, crosslen);
                string idv2part = GeneGroup[gid][gb2num].geneList.substr(0, crosslen);

                GeneGroup[gid][i].geneList.replace(0, crosslen, idv2part);
                GeneGroup[gid][gb2num].geneList.replace(0, crosslen, idv1part);
            }
        }
    }
}

/*
 * 按照变异概率Pm对种群中的个体进行变异
 */
void GeneSFlow::IdvVariation(int gid)
{
    for(int i = 0; i < M; i++)
    {
        double tmpnum = GSFRandom(1000)/1000.0;
        if(tmpnum <= Pm)
        {
            int varLoc = GSFRandom(genelen);
            char c = GeneGroup[gid][i].geneList[varLoc];
            c = '0' + (c - '0' + 1) % 2;
            GeneGroup[gid][i].geneList[varLoc] = c;
        }
    }
}

/*
 * 按比例选择个体到种群1
 */
void GeneSFlow::ChooseIdvToNext(int cgid, int ngid)
{
    GeneGroup[ngid].clear();

    for(int i = 0; i < M; i++)
    {
        double tmpnum = GSFRandom(1000)/1000.0;
        for(int j = 0; j < M; j++)
        {
            tmpnum -= GeneGroup[cgid][j].prob;
            if(tmpnum < 0)
            {
                GeneGroup[ngid].push_back(GeneGroup[cgid][j]);
                break;
            }
        }
    }
}

/*
 * 计算种群中每个个体的适应度，并按适应度值降序排列
 */
void GeneSFlow::CalIdvFitness(int gid)
{
    double resSum = 0;
    for(int i = 0; i < M; i++)
    {
        double genevalue = GSFDeCode(GeneGroup[gid][i].geneList);
        GeneGroup[gid][i].prob = Fx(genevalue);
        resSum += GeneGroup[gid][i].prob;
    }

    for(int i = 0; i < M; i++)
    {
        GeneGroup[gid][i].prob /= resSum;
    }

    sort(GeneGroup[gid].begin(), GeneGroup[gid].end());
}

/*
 * 解码基因
 */
double GeneSFlow::GSFDeCode(string genelist)
{
    double value = 0;
    for(int i = genelen - 1; i >= 0; i--)
    {
        if(genelist[i] == '1')
        {
            value += pow(2.0, genelen-1-i);
        }
    }

    double res = value * groupRange.length / genecount + groupRange.startValue;
    return res;
}

/*
 * 生成目标长度的个体
 */
string GeneSFlow::GetGeneWithLen(int length)
{
    string res = "";
    char tmp;
    for(int i = 0; i < length; i++)
    {
        tmp = '0' + GSFRandom(2);
        res += tmp;
    }
    return res;
}

/*
 * 产生范围内随机数
 */
int GeneSFlow::GSFRandom(int range)
{
    return rand() % range;
}

/*
 * 打印输出指定种群中的所有个体
 */
void GeneSFlow::PrintGroup(int gid)
{
    for(int i = 0; i < GeneGroup[gid].size(); i++)
    {
        GSFLog(GeneGroup[gid][i].geneList + " ");
        if(i % 5 == 4)
        {
            GSFLog("\n");
        }
    }
}

/*
 * Log方法
 */
void GeneSFlow::GSFLog(string log)
{
    //cout << log << endl;
}

void GeneSFlow::SetProbForCross(double d)
{
    Pc = d;
}

void GeneSFlow::SetProbForMtp(double d)
{
    Pm = d;
}