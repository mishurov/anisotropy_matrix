/* ************************************************************************
 * Copyright 2017 Alexander Mishurov
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


#include <SOP/SOP_Node.h>
#include <GU/GU_Detail.h>
#include <GEO/GEO_PointTree.h>


namespace HDK_AMPlugins {


class SOP_AnisotropyMatrix : public SOP_Node
{
public:
	SOP_AnisotropyMatrix(OP_Network *net, const char *name,
						OP_Operator *op);
	virtual ~SOP_AnisotropyMatrix();

	static OP_Node *myConstructor(OP_Network *net, const char *name,
						OP_Operator *op);
	THREADED_METHOD2(
			SOP_AnisotropyMatrix,
			p_pts > 100,
			compute,
			GU_Detail&, gdp_dublicate,
			GEO_PointTreeGAOffset&, tree)
	void computePartial(
		GU_Detail& gdp_dublicate,
		GEO_PointTreeGAOffset& tree,
		const UT_JobInfo &info);

public:
	static PRM_Template myTemplateList[];
	GA_Size  p_pts;

protected:
	virtual OP_ERROR cookMySop(OP_Context &);
private:
	fpreal smoothing_kernel_radius_;
	fpreal search_radius_;
	fpreal scaling_factor_;
	int particles_threshold_;
	GU_Detail* particles_gdp_;
	GU_Detail* geometry_gdp_;
};


} // HDK_AMPlugins namespace

#endif  // HDKPLUGIN_SOPANISOTROPICMATRIX_H_
