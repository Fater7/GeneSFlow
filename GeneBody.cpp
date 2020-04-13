/*
 * @Author: Fater7 
 * @Date: 2019-07-04 17:07:55 
 * @Last Modified by: Fater7
 * @Last Modified time: 2019-07-04 20:56:28
 */

#include "GeneBody.h"

GeneBody::GeneBody(string tG)
{
    geneList = tG;
    prob = 0.0;
}

bool GeneBody::operator<(const GeneBody& other) const
{
    return prob > other.prob;
}