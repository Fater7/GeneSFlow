/*
 * @Author: Fater7 
 * @Date: 2019-07-04 17:05:37 
 * @Last Modified by: Fater7
 * @Last Modified time: 2019-07-04 20:56:31
 */

#pragma once

#include <string>

using namespace std;

typedef struct GeneBody
{
    GeneBody(string tG);
    string geneList;
    double prob;

    bool operator<(const GeneBody& other) const;
}GeneBody;