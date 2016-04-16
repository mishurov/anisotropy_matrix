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

#include "SOP_anisotropy_matrix.h"

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
                                     0));
}

static PRM_Name krnl_name("kernel", "Kernel Radius");
static PRM_Name srch_name("search", "Search Radius");
static PRM_Name scl_name("scale", "Scale Addition");
static PRM_Name thrs_name("threshold", "Particles Threshold");

static PRM_Default krnl_default(1);
static PRM_Default srch_default(2);
static PRM_Default scl_default(0.1);
static PRM_Default thrs_default(6);

static PRM_Range unsigned_range(PRM_RANGE_RESTRICTED, 0);
static PRM_Range unsigned_10_range(PRM_RANGE_UI, 0, PRM_RANGE_RESTRICTED, 10);
static PRM_Range unsigned_50_range(PRM_RANGE_UI, 0, PRM_RANGE_RESTRICTED, 50);

PRM_Template
SOP_AnisotropyMatrix::myTemplateList[] = {
  PRM_Template(PRM_FLT,  1, &krnl_name, &krnl_default, 0, &unsigned_10_range),
  PRM_Template(PRM_FLT,  1, &srch_name, &srch_default, 0, &unsigned_10_range),
  PRM_Template(PRM_FLT,  1, &scl_name, &scl_default, 0, &unsigned_range),
  PRM_Template(PRM_INT,  1, &thrs_name, &thrs_default, 0, &unsigned_50_range),
  PRM_Template(),
};

OP_Node *
SOP_AnisotropyMatrix::myConstructor(OP_Network *net,
                                    const char *name,
                                    OP_Operator *op)
{
    return new SOP_AnisotropyMatrix(net, name, op);
}

SOP_AnisotropyMatrix::SOP_AnisotropyMatrix(OP_Network *net,
                                           const char *name,
                                           OP_Operator *op)
: SOP_Node(net, name, op) {}

SOP_AnisotropyMatrix::~SOP_AnisotropyMatrix() {}

OP_ERROR
SOP_AnisotropyMatrix::cookMySop(OP_Context &context)
{
  if (lockInputs(context) >= UT_ERROR_ABORT)
    return error();

  setupLocalVars();

  if (error() < UT_ERROR_ABORT) {
    UT_AutoInterrupt progress("Calculating matrices");

    GU_Detail *particles_gdp = (GU_Detail *)inputGeo(0, context);
    GU_Detail *geometry_gdp = (GU_Detail *)inputGeo(1, context);

    gdp->clearAndDestroy();
    GU_Detail gdp_dublicate(true);
    gdp_dublicate.copy(*geometry_gdp,
                       GEO_COPY_ONCE,
                       false,
                       true,
                       GA_DATA_ID_BUMP);

    // MATRIX ATTRIBUTES
    fpreal smoothing_kernel_radius = evalFloat("kernel", 0, 0);
    fpreal search_radius = evalFloat("search", 0, 0);
    fpreal scale_addition = evalFloat("scale", 0, 0);
    unsigned particles_threshold = evalInt("threshold", 0, 0);

    GEO_PointTreeGAOffset tree;
    tree.build(particles_gdp, NULL);

    GA_Size p_pts = particles_gdp->getNumPoints();

    for (unsigned ii=0; ii<p_pts; ii++) {
      UT_Vector3 particle_pos = particles_gdp->getPos3(ii);

      // Close particles indices
      GEO_PointTreeGAOffset::IdxArrayType close_particles_indices;
      tree.findAllCloseIdx(particle_pos, 
                           search_radius,
                           close_particles_indices);

      // NOTE from HDK 13:
      // entries() will be renamed to size() in a future version
      unsigned close_particles_count = close_particles_indices.entries();

      UT_Matrix3 anisotropy_matrix;

      if (close_particles_count > 0) {

        // Calculation of weighted mean
        UT_Vector3 weighted_mean(0, 0, 0);
        fpreal weight = 0;
        fpreal weighting_function = 0;
        UT_Vector3 weighted_position(0, 0, 0);

        for (unsigned i = 0; i < close_particles_count; i++) {
          UT_Vector3 close_particle_pos = 
          particles_gdp->getPos3(close_particles_indices(i));

          UT_Vector3 distance = particle_pos - close_particle_pos;
          weight = 1 - std::pow((distance.length() / search_radius), 3);

          weighting_function += weight;
          weighted_position += close_particle_pos * weight;
        }

        if (weighting_function != 0)
          weighted_mean = weighted_position/weighting_function;
        else
          weighted_mean = weighted_position;

        // Calculation of covariance matrix and SVD -- example code provided by ndickson, thank you!
        UT_Matrix3D covariance_matrix(0); 
        for (unsigned i = 0; i < close_particles_count; i++) { 
            UT_Vector3 close_particle_pos = 
                particles_gdp->getPos3(close_particles_indices(i)); 
            UT_Vector3 weighted_distance = close_particle_pos - weighted_mean; 
            UT_Vector3 distance = particle_pos - close_particle_pos; 
            weight = 1 - std::pow((distance.length() / search_radius), 3); 
            UT_Vector3 weighted = weight * weighted_distance; 
            // Only 6 unique components, since symmetric 
            covariance_matrix(0,0) += weighted(0)*weighted_distance(0); 
            covariance_matrix(0,1) += weighted(0)*weighted_distance(1); 
            covariance_matrix(0,2) += weighted(0)*weighted_distance(2); 
            covariance_matrix(1,1) += weighted(1)*weighted_distance(1); 
            covariance_matrix(1,2) += weighted(1)*weighted_distance(2); 
            covariance_matrix(2,2) += weighted(2)*weighted_distance(2); 
        } 

        // Copy symmetric components 
        covariance_matrix(1,0) = covariance_matrix(0,1); 
        covariance_matrix(2,0) = covariance_matrix(0,2); 
        covariance_matrix(2,1) = covariance_matrix(1,2); 

        if (weighting_function != 0) 
            covariance_matrix /= weighting_function; 

        UT_Matrix3D rotation_matrix; 
        UT_Matrix3D diagonal_matrix; // diagonal will hold eigenvalues 

        // This is probably overkill; there are likely faster algorithms, but this should work, up to some tolerance. 
        covariance_matrix.diagonalizeSymmetric(rotation_matrix,diagonal_matrix); 

        UT_Vector3D eigen_values_vector(diagonal_matrix(0,0), diagonal_matrix(1,1), diagonal_matrix(2,2));
        // End SVD/CovarianceMatrix

        // Particles threshold
        if (close_particles_count > particles_threshold) {
          for (unsigned i = 0; i < 3; i++)
            eigen_values_vector(i) = 2 * eigen_values_vector(i) +
                                     scale_addition;
        } else {
          eigen_values_vector(0) =
          eigen_values_vector(1) =
          eigen_values_vector(2) = 1;
        }

        // Convert to HDK matrices
        UT_Matrix3 rotation_matrix_hdk(rotation_matrix(0, 0),
                                       rotation_matrix(0, 1),
                                       rotation_matrix(0, 2),
                                       rotation_matrix(1, 0),
                                       rotation_matrix(1, 1),
                                       rotation_matrix(1, 2),
                                       rotation_matrix(2, 0),
                                       rotation_matrix(2, 1),
                                       rotation_matrix(2, 2));

        UT_Matrix3 eigen_values_matrix_hdk(eigen_values_vector(0), 0, 0,
                                           0, eigen_values_vector(1), 0,
                                           0, 0, eigen_values_vector(2));

        rotation_matrix_hdk.invert();
        UT_Matrix3 rotation_matrix_transpose_hdk = rotation_matrix_hdk;
        rotation_matrix_transpose_hdk.transpose();

        // Compose anisotropy matrix
        anisotropy_matrix = rotation_matrix_hdk *
                            eigen_values_matrix_hdk *
                            rotation_matrix_transpose_hdk;
        anisotropy_matrix *= 1 / smoothing_kernel_radius;
      }

      // Calculate new transormations for dublicate geometry
      GA_Size g_pts = geometry_gdp->getNumPoints();


    for (unsigned ee=0; ee < g_pts; ee++) {
      UT_Vector3 geometry_pos = geometry_gdp->getPos3(ee);
      if (close_particles_count > 0)
        geometry_pos.colVecMult(anisotropy_matrix);

      gdp_dublicate.setPos3(ee, geometry_pos + particle_pos);
    }

    // Add geometry copy to final geometry
    gdp->copy(gdp_dublicate,
              GEO_COPY_ADD,
              true,
              true,
              GA_DATA_ID_BUMP);

    }
  }

  unlockInputs();
  resetLocalVarRefs();

  return error();
}
