/* ************************************************************************
 * Copyright 2013 Alexander Mishurov
 * 
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 * 
 * http://www.apache.org/licenses/LICENSE-2.0
 * 
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 * ************************************************************************/

#ifndef HDKPLUGIN_SOPANISOTROPICMATRIX_H_
#define HDKPLUGIN_SOPANISOTROPICMATRIX_H_

#include <omp.h>
#include <SYS/SYS_Math.h>
#include <UT/UT_DSOVersion.h>
#include <UT/UT_Interrupt.h>
#include <UT/UT_Matrix4.h>
#include <OP/OP_Operator.h>
#include <OP/OP_OperatorTable.h>
#include <SOP/SOP_Node.h>
#include <PRM/PRM_Include.h>
#include <GU/GU_Detail.h>
#include <GEO/GEO_PointTree.h>

namespace HDK_AMPlugins {

class SOP_AnisotropyMatrix : public SOP_Node
{
public:
  SOP_AnisotropyMatrix(OP_Network *, const char *, OP_Operator *);
  virtual ~SOP_AnisotropyMatrix();
  static OP_Node *myConstructor(OP_Network *, const char *, OP_Operator *);
public:
  static PRM_Template myTemplateList[];

protected:
  virtual OP_ERROR cookMySop(OP_Context &);
};

} // HDK_AMPlugins namespace

#endif  // HDKPLUGIN_SOPANISOTROPICMATRIX_H_
