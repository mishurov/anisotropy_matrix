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


#include "SOP_anisotropy_matrix.h"


#include <GA/GA_Iterator.h>
#include <OP/OP_AutoLockInputs.h>
#include <OP/OP_Operator.h>
#include <OP/OP_OperatorTable.h>
#include <PRM/PRM_Include.h>
#include <CH/CH_LocalVariable.h>
#include <UT/UT_DSOVersion.h>
#include <UT/UT_Interrupt.h>
#include <UT/UT_Matrix3.h>
#include <UT/UT_Vector3.h>


using namespace HDK_AMPlugins;


void
newSopOperator(OP_OperatorTable *table)
{
	table->addOperator(new OP_Operator("hdk_anisotropy_matrix",
		"AnisotropyMatrix",
		SOP_AnisotropyMatrix::myConstructor,
		SOP_AnisotropyMatrix::myTemplateList,
		2,
		2,
		NULL));
}


static PRM_Name kernel_name("kernel", "Kernel Radius");
static PRM_Default kernel_default(1);

static PRM_Name search_name("search", "Search Radius");
static PRM_Default search_default(2);

static PRM_Name scale_name("scale", "Scaling Factor");
static PRM_Default scale_default(1);

static PRM_Name threshold_name("threshold", "Particles Threshold");
static PRM_Default threshold_default(6);


static PRM_Range range10(PRM_RANGE_UI, 0, PRM_RANGE_RESTRICTED, 10);
static PRM_Range range50(PRM_RANGE_UI, 0, PRM_RANGE_RESTRICTED, 50);


PRM_Template
SOP_AnisotropyMatrix::myTemplateList[] = {
	PRM_Template(PRM_FLT, 1, &kernel_name, &kernel_default, 0, &range10),
	PRM_Template(PRM_FLT, 1, &search_name, &search_default, 0, &range10),
	PRM_Template(PRM_FLT, 1, &scale_name, &scale_default, 0, &range10),
	PRM_Template(
		PRM_INT, 1, &threshold_name, &threshold_default, 0, &range50
	),
	PRM_Template(),
};


OP_Node *
SOP_AnisotropyMatrix::myConstructor(OP_Network *net, const char *name,
							OP_Operator *op)
{
	return new SOP_AnisotropyMatrix(net, name, op);
}


SOP_AnisotropyMatrix::SOP_AnisotropyMatrix(OP_Network *net, const char *name,
							OP_Operator *op)
	: SOP_Node(net, name, op) {}


SOP_AnisotropyMatrix::~SOP_AnisotropyMatrix() {}


void
SOP_AnisotropyMatrix::computePartial(GU_Detail& gdp_dublicate,
			GEO_PointTreeGAOffset& tree, const UT_JobInfo &info)
{
int i, n;
for (info.divideWork(p_pts, i, n); i < n; i++) {
	UT_Vector3 particle_pos = particles_gdp_->getPos3(i);

	// Close particles indices
	GEO_PointTreeGAOffset::IdxArrayType close_particles_indices;
	tree.findAllCloseIdx(
		particle_pos, search_radius_, close_particles_indices
	);

	int close_particles_count = close_particles_indices.size();
	
	UT_Matrix3 anisotropy_matrix;

	if (close_particles_count > 0) {
		// Compute a weighted mean
		UT_Vector3 weighted_mean(0, 0, 0);
		fpreal weight = 0;
		fpreal weighting_function = 0;
		UT_Vector3 weighted_position(0, 0, 0);

		for (unsigned j = 0; j < close_particles_count; j++) {
			UT_Vector3 close_particle_pos = 
				particles_gdp_->getPos3(
					close_particles_indices(j)
				);

			UT_Vector3 distance = particle_pos - close_particle_pos;
			weight = 1 - pow(
				(distance.length() / search_radius_), 3
			);

			weighting_function += weight;
			weighted_position += close_particle_pos * weight;
		}

		weighted_mean = weighting_function != 0 ?
			weighted_position / weighting_function :
				weighted_position;

		// Compute a covariance matrix
		UT_MatrixT<float> weighted_distance_column;
		weighted_distance_column.resize(3, 1);

		UT_MatrixT<float> weighted_distance_row;
		weighted_distance_row.resize(1, 3);

		UT_MatrixT<float> covariance_matrix;
		covariance_matrix.resize(3, 3);
		covariance_matrix.zero();

		for (unsigned j = 0; j < close_particles_count; j++) {
			UT_Vector3 close_particle_pos =
				particles_gdp_->getPos3(
					close_particles_indices(j)
				);
			UT_Vector3 weighted_distance =
					close_particle_pos - weighted_mean;
			
			for (int k = 0; k < 3; k++) {
				weighted_distance_column(k, 0) =
							weighted_distance(k);
				weighted_distance_row(0, k) =
							weighted_distance(k);
			}

			UT_Vector3 distance = particle_pos - close_particle_pos;
			weight = 1 - pow(
				(distance.length() / search_radius_), 3
			);

			UT_MatrixT<float> res;
			res.resize(3, 3);
			weighted_distance_column.postMult(
				weighted_distance_row, res
			);
			covariance_matrix.addScaledMatrix(res, weight);
		}

		if (weighting_function != 0) {
			covariance_matrix.setAndScale(
				covariance_matrix, 1 / weighting_function
			);
		}

		// Singular Value Decomposition to get matrices
		UT_Matrix3 covariance_matrix3;
		covariance_matrix.getSubmatrix3(covariance_matrix3, 0, 0);

		UT_Matrix3 rotation_matrix;
		UT_Matrix3 eigen_values_matrix;
		UT_Matrix3 rotation_matrix_transpose;
		covariance_matrix3.svdDecomposition(
			rotation_matrix,
			eigen_values_matrix,
			rotation_matrix_transpose
		);

		rotation_matrix_transpose = rotation_matrix;
		rotation_matrix_transpose.transpose();
		// Particles threshold
		if (close_particles_count > particles_threshold_) {
			eigen_values_matrix *= scaling_factor_;
		} else {
			eigen_values_matrix.identity();
		}

		// Compose anisotropy matrix
		anisotropy_matrix = rotation_matrix *
					eigen_values_matrix *
						rotation_matrix_transpose;
		anisotropy_matrix *= smoothing_kernel_radius_;
	}

	// Calculate new transormations for dublicate geometry
	GA_Size g_pts = geometry_gdp_->getNumPoints();

	// Shared between threads data
	{
		UT_AutoJobInfoLock a(info);

		for (unsigned j = 0; j < g_pts; j++) {
			UT_Vector3 geometry_pos = geometry_gdp_->getPos3(j);
			if (close_particles_count > 0) {
				geometry_pos.colVecMult(anisotropy_matrix);
			}

			gdp_dublicate.setPos3(j, geometry_pos + particle_pos);
		}

		// Add geometry copy to the final geometry
		gdp->copy(gdp_dublicate,
			  GEO_COPY_ADD,
			  true,
			  true,
			  GA_DATA_ID_BUMP);
	}

} // endfor
}


OP_ERROR
SOP_AnisotropyMatrix::cookMySop(OP_Context &context)
{
	OP_AutoLockInputs inputs(this);
	if (inputs.lock(context) >= UT_ERROR_ABORT)
		return error();

	setupLocalVars();

	if (error() < UT_ERROR_ABORT) {
		UT_AutoInterrupt progress("Calculating matrices");

		particles_gdp_ = (GU_Detail *)inputGeo(0, context);
		geometry_gdp_ = (GU_Detail *)inputGeo(1, context);

		gdp->clearAndDestroy();
		GU_Detail gdp_dublicate(true);
		gdp_dublicate.copy(
			*geometry_gdp_,
			GEO_COPY_ONCE,
			false,
			true,
			GA_DATA_ID_BUMP
		);

		// attributes
		smoothing_kernel_radius_ = evalFloat("kernel", 0, 0);
		search_radius_ = evalFloat("search", 0, 0);
		scaling_factor_ = evalFloat("scale", 0, 0);
		particles_threshold_ = evalInt("threshold", 0, 0);

		GEO_PointTreeGAOffset tree;
		tree.build(particles_gdp_, NULL);

		p_pts = particles_gdp_->getNumPoints();
	
		// multithreaded evaluation
		compute(gdp_dublicate, tree);

		// single threaded
		//computeNoThread(gdp_dublicate, tree);
	}

	resetLocalVarRefs();

	return error();
}
