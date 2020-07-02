#include "diffgeom.h"

//=====================  BODY  ==========================

diffgeom::diffgeom(){}

diffgeom::diffgeom(string input){
  RunMode = input;
}

diffgeom::~diffgeom(){}

void diffgeom::get_coordinates(vector<vector<double> > comp_coords, vector<vector<double> > phys_coords){
  comp_nodes.resize(Nx*Ny*Nz);
  for(unsigned int i=0; i<comp_nodes.size(); i++){
    comp_nodes[i].resize(Dim);
  }
  phys_nodes.resize(Nx*Ny*Nz);
  for(unsigned int i=0; i<phys_nodes.size(); i++){
    phys_nodes[i].resize(Dim);
  }
  if(comp_nodes.size()==comp_coords.size() && phys_nodes.size()==phys_coords.size() && comp_nodes[0].size()==phys_nodes[0].size() && comp_nodes[0].size()==comp_coords[0].size() && phys_nodes[0].size()==phys_coords[0].size() && comp_nodes[0].size()==phys_nodes[0].size()){
    for(unsigned int i=0; i<comp_nodes.size(); i++){
      for(unsigned int j=0; j<comp_nodes[i].size(); j++){
        comp_nodes[i][j] = comp_coords[i][j]; 
      }
    }
    for(unsigned int i=0; i<comp_nodes.size(); i++){
      for(unsigned int j=0; j<comp_nodes[i].size(); j++){
        phys_nodes[i][j] = phys_coords[i][j]; 
      }
    }
  }
  else {
    cout << "Error: mismatch in vector dimensions. Control the assignments you made.";
    exit(1);
  }
}

void diffgeom::compute_covbase(){
  
  cov_base.resize(Nx*Ny*Nz);
  for(unsigned int i=0; i<cov_base.size(); i++){
    cov_base[i].resize(Dim*Dim);
  }
  
  // dev_order is the derivative order, 1 for the current case
  // dev_order is the equivalent derivative order, needed to compute the stencil dimension given the accuracy
  int dev_order, dev_order_eq;
  dev_order = 1;
  if(dev_order%2==0 && dev_order!=0){
      dev_order_eq = dev_order-1;
  }else{
      dev_order_eq = dev_order;
  }
  // Stencil dimension
  int n_stencil;
  if(dev_order!=0){
        n_stencil = dev_order_eq + A;
      }else{
        n_stencil = 1;
      }
  int transl = (int)(n_stencil-1)/2;                                 // Index translation to fit stencil into grid
  
  // Initialize the computation of the covariant vectors
  int start_index;                                                   // Index of first stencil's node
  int final_index;                                                   // Index of last stencil's node
  int curr_n_stencil;                                                // Current stencil dimension
  vector<int> curr_indexes;                                          // Current indexes of grid points [i,j,k]
  vector<int> grid_dim;                                              // Number of grid points in each direction [Nx, Ny, Nz]
  grid_dim.resize(3);
  grid_dim[0] = Nx;
  grid_dim[1] = Ny;
  grid_dim[2] = Nz;
  vector<vector<int> > stencil;                                      // Grid indexes of stencil nodes [[i_(1),j_(1),k_(1)], ..., [i_(n_stencil),j_(n_stencil),k_(n_stencil)]]
  vector<double> weights;                                            // Weights of finite difference approximation
  vector<double> xs;                                                 // Stencil reference coordinates
  vector<double> fs;                                                 // Function values at stencil grid points
  double z;                                                          // Derivative's evaluation point
  
  // Compute the covariant vectors
  //--- Loop starts ---
  for(int cov=0; cov<Dim; cov++){                                    // Loop through the vector (2 or 3, depending on space Dim)
    for(int comp=0; comp<Dim; comp++){                               // Loop through each covariant vector's components (2 or 3, depending on space Dim)
      // Loop through the grid
      for(int k=0; k<Nz; k++){
        for(int j=0; j<Ny; j++){
          for(int i=0; i<Nx; i++){
            
            // Re-set to 0, to avoid errors
            start_index = 0;
            final_index = 0;
            curr_n_stencil = 0;
            curr_indexes.clear();
            stencil.clear();
            weights.clear();
            xs.clear();
            fs.clear();
            z=0;
            
            curr_n_stencil = n_stencil;
            //index = i + j*Nx + k*Nx*Ny;                               // Compute the helical index. Note: in 2D, Nz = 0; thus, k =[0] for which k*Nx*Ny = 0
            curr_indexes.resize(3);
            curr_indexes[0] = i;
            curr_indexes[1] = j;
            curr_indexes[2] = k;
            start_index = curr_indexes[cov] - transl;
            final_index = curr_indexes[cov] + transl;
            if(start_index<0 && final_index<grid_dim[cov]){                // Leftmost stencil node outside index range (smaller than 0)
              start_index = 0;
            }
            else if(start_index>=0 && final_index>=grid_dim[cov]){         // Rightmost stencil node outside index range (greater or equal than grid size)
              start_index -= (final_index-(grid_dim[cov]-1));
            }
            else if(start_index<0 && final_index>=grid_dim[cov]){         // Leftmost and rightmost stencil node outside index range ==> reduce stencil size (accuracy is reduced)
              curr_n_stencil = grid_dim[cov];
              start_index = 0;
            }
            stencil.resize(curr_n_stencil);
            for(int n=0; n<stencil.size(); n++){
              stencil[n].resize(3);
              for(int m=0; m<stencil[n].size(); m++){
                stencil[n][m] = curr_indexes[m];
              }
            }
            for(int n=0; n<stencil.size(); n++){
              stencil[n][cov] = start_index + n;
            }
            xs.resize(stencil.size());
            for(int n=0; n<stencil.size(); n++){
              xs[n] = comp_nodes[stencil[n][0] + stencil[n][1]*Nx + stencil[n][2]*Nx*Ny][cov];
            }
            fs.resize(stencil.size());
            for(int n=0; n<stencil.size(); n++){
              fs[n] = phys_nodes[stencil[n][0] + stencil[n][1]*Nx + stencil[n][2]*Nx*Ny][comp];
            }
            z = comp_nodes[i + j*Nx + k*Nx*Ny][cov];
            
            finitediff differences(dev_order,z,xs,fs);
            differences.compute_weights();
            differences.compute_derivative();
            cov_base[i + j*Nx + k*Nx*Ny][comp + cov*Dim] = differences.provide_derivative();
            
          }
        }
      }
    }
  }
  //--- Loop finished ---
}

void diffgeom::compute_covmetric(){
  cov_metric.resize(Nx*Ny*Nz);
  for(unsigned int i=0; i<cov_metric.size(); i++){
    cov_metric[i].resize(Dim*Dim);
  }
  for(int k=0; k<Nz; k++){
    for(int j=0; j<Ny; j++){
      for(int i=0; i<Nx; i++){
        for(int left_cov=0; left_cov<Dim; left_cov++){                                 
          for(int right_cov=0; right_cov<Dim; right_cov++){
            cov_metric[i + j*Nx + k*Nx*Ny][right_cov + left_cov*Dim] = 0;
            for(int n=0; n<Dim; n++){
              cov_metric[i + j*Nx + k*Nx*Ny][right_cov + left_cov*Dim] += cov_base[i + j*Nx + k*Nx*Ny][left_cov*Dim+n]*cov_base[i + j*Nx + k*Nx*Ny][right_cov*Dim+n];
            }
          }
        }
      }
    }
  }
}

vector<double> diffgeom::compute_metricdet(){
  vector<double> metricdet;
  metricdet.resize(Nx*Ny*Nz);
  
  vector<vector<double> > A;
  
  for(int k=0; k<Nz; k++){
    for(int j=0; j<Ny; j++){
      for(int i=0; i<Nx; i++){
        A.resize(Dim);
        for(int n=0; n<Dim; n++){
          A[n].resize(Dim)
          for(int m=0; m<Dim; m++){
            A[n][m] = cov_metric[i + j*Nx + k*Nx*Ny][m + n*Dim];
          }
        }
        matrix mat(A);
        mat.det();
        metricdet[i + j*Nx + k*Nx*Ny] = mat.provide_det();
        A.clear();
      }
    }
  }
  
  return metricdet;
}

void diffgeom::compute_contrbase(){
  contr_base.resize(Nx*Ny*Nz);
  for(unsigned int i=0; i<contr_base.size(); i++){
    contr_base[i].resize(Dim*Dim);
  }
}

void diffgeom::compute_contrmetric(){
  contr_metric.resize(Nx*Ny*Nz);
  for(unsigned int i=0; i<contr_metric.size(); i++){
    contr_metric[i].resize(Dim*Dim);
  }
  for(int k=0; k<Nz; k++){
    for(int j=0; j<Ny; j++){
      for(int i=0; i<Nx; i++){
        for(int left_contr=0; left_contr<Dim; left_contr++){                                 
          for(int right_contr=0; right_contr<Dim; right_contr++){
            contr_metric[i + j*Nx + k*Nx*Ny][right_contr + left_contr*Dim] = 0;
            for(int n=0; n<Dim; n++){
              contr_metric[i + j*Nx + k*Nx*Ny][right_contr + left_contr*Dim] += contr_base[i + j*Nx + k*Nx*Ny][left_contr*Dim+n]*contr_base[i + j*Nx + k*Nx*Ny][right_contr*Dim+n];
            }
          }
        }
      }
    }
  }
}

void diffgeom::compute_firstChristoffel(){
  first_christoffel.resize(Nx*Ny*Nz);
  for(unsigned int i=0; i<first_christoffel.size(); i++){
    first_christoffel[i].resize(Dim*Dim*Dim);
  }
}

void diffgeom::compute_secondChristoffel(){
  second_christoffel.resize(Nx*Ny*Nz);
  for(unsigned int i=0; i<second_christoffel.size(); i++){
    second_christoffel[i].resize(Dim*Dim*Dim);
  }
}

void diffgeom::compute_contrRiemann(){
  contr_Riemann.resize(Nx*Ny*Nz);
  for(unsigned int i=0; i<contr_Riemann.size(); i++){
    contr_Riemann[i].resize(Dim*Dim*Dim*Dim);
  }
}

void diffgeom::compute_covRiemann(){
  cov_Riemann.resize(Nx*Ny*Nz);
  for(unsigned int i=0; i<cov_Riemann.size(); i++){
    cov_Riemann[i].resize(Dim*Dim*Dim*Dim);
  }
}

void diffgeom::compute_covRicci(){
  cov_Ricci.resize(Nx*Ny*Nz);
  for(unsigned int i=0; i<cov_Ricci.size(); i++){
    cov_Ricci[i].resize(Dim*Dim);
  }
}

void diffgeom::compute_Ricciscalar(){
  Ricci_scalar.resize(Nx*Ny*Nz);

}

vector<vector<double> > diffgeom::compute_invmetricdet(){
  vector<vector<double> > invmetricdet;
  invmetricdet.resize(Nx*Ny*Nz);
  
  return invmetricdet;
}

vector<vector<double> > diffgeom::compute_covveccomp(vector<vector<double> >){
  vector<vector<double> > covveccomp;
  covveccomp.resize(Nx*Ny*Nz);
  for(unsigned int i=0; i<covveccomp.size(); i++){
    covveccomp[i].resize(Dim*Dim);
  }
  
  return covveccomp;
}

vector<vector<double> > diffgeom::compute_contrveccomp(vector<vector<double> >){
  vector<vector<double> > contrveccomp;
  contrveccomp.resize(Nx*Ny*Nz);
  for(unsigned int i=0; i<contrveccomp.size(); i++){
    contrveccomp[i].resize(Dim*Dim);
  }
  
  return contrveccomp;
}

vector<vector<double> > diffgeom::compute_vecfromcovcomp(vector<vector<double> >){
  vector<vector<double> > vec;
  vec.resize(Nx*Ny*Nz);
  for(unsigned int i=0; i<vec.size(); i++){
    vec[i].resize(Dim*Dim);
  }
  
  return vec;
}

vector<vector<double> > diffgeom::compute_vecfromcontrcomp(vector<vector<double> >){
  vector<vector<double> > vec;
  vec.resize(Nx*Ny*Nz);
  for(unsigned int i=0; i<vec.size(); i++){
    vec[i].resize(Dim*Dim);
  }
  
  return vec;
}

vector<vector<double> > diffgeom::provide_compdomain(){
  return comp_nodes;
}

vector<vector<double> > diffgeom::provide_physdomain(){
  return phys_nodes;
}

vector<vector<double> > diffgeom::provide_covbase(){
  return cov_base;
}

vector<vector<double> > diffgeom::provide_contrbase(){
  return contr_base;
}

vector<vector<double> > diffgeom::provide_covmetric(){
  return cov_metric;
}

vector<vector<double> > diffgeom::provide_contrmetric(){
  return contr_metric;
}

vector<vector<double> > diffgeom::provide_firstChristoffel(){
  return first_christoffel;
}

vector<vector<double> > diffgeom::provide_secondChristoffel(){
  return second_christoffel;
}

vector<vector<double> > diffgeom::provide_contrRiemann(){
  return contr_Riemann;
}

vector<vector<double> > diffgeom::provide_covRiemann(){
  return cov_Riemann;
}

vector<vector<double> > diffgeom::provide_covRicci(){
  return cov_Ricci;
}

double diffgeom::provide_Ricciscalar(){
  return Ricci_scalar;
}