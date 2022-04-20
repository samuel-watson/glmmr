#include <RcppArmadillo.h>
#include "glmmrmath.h"
#include "glmmrmatrix.h"
using namespace Rcpp;
using namespace arma;

#ifdef _OPENMP
#include <omp.h>
#endif

// [[Rcpp::depends(RcppArmadillo)]]


class HillClimbing {
private:
  const arma::field<arma::vec> C_list_; 
  arma::field<arma::mat> X_list_; 
  arma::field<arma::mat> X_list_sub_; 
  const arma::field<arma::mat> D_list_;
  const arma::field<arma::mat> X_all_list_; 
  const arma::field<arma::mat> Z_all_list_; 
  const arma::mat W_all_diag_;
  
  arma::uvec all_idx_;
  const arma::vec weights_; // weights for each design
  arma::vec count_in_design_;
  arma::uvec max_obs_;
  arma::uvec curr_obs_;
  const arma::uword nlist_; // number of designs
  arma::uword any_fix_;
  const arma::uword n_; //size of the design to find
  arma::uword k_;
  arma::umat track_swaps_;
  
public:
  arma::uvec idx_in_; 
  arma::uvec idx_in_rm_;
  arma::uvec idx_in_sub_;
  //arma::uvec idx_out_; 
  arma::vec new_val_vec_; //value of objective function for each separate design
  arma::vec best_val_vec_; //value of objective function best
  double val_; // overall best value
  double new_val_; // new value
  double rm_val_;
  arma::uword s_;
  double r_;
  double b_;
  
  arma::cube A_list_; // inverse sigma matrices
  arma::cube A_list_sub_; // holder of 
  arma::cube rm1A_list_; // inverse sigma matrices with one removed - initialised to minus one but now needs to resize
  arma::field<arma::mat> M_list_; // field of information matrices
  arma::field<arma::mat> M_list_sub_; // field of information matrices
  const arma::uvec nfix_; //the indexes of the experimental conditions to keep
  const arma::uword rd_mode_; // robust designs mode: 1 == weighted, 2 == minimax.
  
  bool trace_;
  
public:
  HillClimbing(arma::uvec idx_in, 
               arma::field<arma::vec> C_list, 
               arma::field<arma::mat> X_list, 
               arma::field<arma::mat> Z_list, 
               arma::field<arma::mat> D_list,
               arma::mat w_diag,
               arma::uvec max_obs,
               arma::vec weights,
               arma::uword any_fix,
               arma::uvec nfix, 
               arma::uword s=10,
               double r = 0.01,
               double b = 1,
               arma::uword rd_mode = 0, 
               bool trace=false) :
  C_list_(C_list), 
  X_list_(weights.n_elem),
  X_list_sub_(weights.n_elem),
  D_list_(D_list),
  X_all_list_(X_list),
  Z_all_list_(Z_list),
  W_all_diag_(w_diag),
  all_idx_(idx_in.n_elem),
  weights_(weights), 
  count_in_design_(X_list(0,0).n_rows,fill::zeros),
  max_obs_(max_obs),
  curr_obs_(X_list(0,0).n_rows,fill::zeros),
  nlist_(weights_.n_elem),
  any_fix_(any_fix),
  n_(idx_in.n_elem), 
  k_(X_list(0,0).n_rows),
  track_swaps_(n_,2,fill::zeros),
  idx_in_(idx_in),
  idx_in_rm_(n_-1,arma::fill::zeros),
  idx_in_sub_(idx_in),
  new_val_vec_(nlist_, arma::fill::zeros), 
  best_val_vec_(nlist_, arma::fill::zeros), 
  val_(0.0), 
  new_val_(0.0),
  rm_val_(0.0),
  s_(s),
  r_(r),
  b_(b),
  A_list_(n_,n_,nlist_,fill::zeros), 
  A_list_sub_(n_,n_,nlist_,fill::zeros),
  rm1A_list_(n_-1,n_-1,nlist_,fill::zeros),
  M_list_(nlist_),
  M_list_sub_(nlist_),
  nfix_(nfix), 
  rd_mode_(rd_mode), 
  trace_(trace) {
    build_XZ();
    Update_M_list(true);
  }
  
  
  arma::uvec join_idx(const arma::uvec &idx, arma::uword elem) {
    return arma::join_cols(idx, arma::uvec({elem}));
  }
  
  void build_XZ(){
    idx_in_ = sort(idx_in_);
    for(arma::uword i=0; i<n_; i++){
      curr_obs_(idx_in_(i))++;
    }
    idx_in_sub_ = sort(idx_in_);
    for(arma::uword j=0; j<nlist_;j++){
      arma::mat X(n_,X_all_list_(j,0).n_cols,fill::zeros);
      arma::mat Z(n_,Z_all_list_(j,0).n_cols,fill::zeros);
      arma::vec w_diag(n_);
      for(arma::uword k=0; k<n_;k++){
        X.row(k) = X_all_list_(j,0).row(idx_in_(k));
        Z.row(k) = Z_all_list_(j,0).row(idx_in_(k));
        w_diag(k) = W_all_diag_(idx_in_(k),j);
        if(j==0)count_in_design_(idx_in_(k)) += 1;
      }
      X_list_(j,0) = X;
      X_list_sub_(j,0) = X;
      arma::vec idx = arma::linspace(0, n_-1, n_);
      all_idx_ = arma::conv_to<arma::uvec>::from(idx);
      arma::mat tmp = Z*D_list_(j,0)*Z.t();
      tmp.diag() += w_diag;
      A_list_.slice(j) = tmp.i();
      A_list_sub_.slice(j) = tmp.i();
    }
  }
  
  double comp_c_obj_fun() {
    for (arma::uword j = 0; j < nlist_; ++j) {
      new_val_vec_(j) = c_obj_fun(arma::conv_to<arma::mat>::from(M_list_sub_(j,0)), 
                   C_list_(j,0));
    }
    return rd_mode_ == 1 ? arma::dot(new_val_vec_, weights_) : new_val_vec_.max();
  }
  
  double direct_deriv(){
    arma::vec res(nlist_);
    for (arma::uword j = 0; j < nlist_; ++j) {
      res(j) = c_d_deriv(arma::conv_to<arma::mat>::from(M_list_(j,0)), 
          arma::conv_to<arma::mat>::from(M_list_sub_(j,0)), 
          C_list_(j,0));
    }
    return rd_mode_ == 1 ? arma::dot(res, weights_) : res.max();
  }
  
  double fnorm(){
    arma::vec res(nlist_);
    for (arma::uword j = 0; j < nlist_; ++j) {
      res(j) = norm(arma::conv_to<arma::mat>::from(M_list_(j,0))- 
        arma::conv_to<arma::mat>::from(M_list_sub_(j,0)),"fro");
    }
    return rd_mode_ == 1 ? arma::dot(res, weights_) : res.max();
  }
  
  double grad_norm(){
    arma::vec res(nlist_);
    for (arma::uword j = 0; j < nlist_; ++j) {
      res(j) =c_deriv(arma::conv_to<arma::mat>::from(M_list_(j,0)),
          C_list_(j,0));
    }
    return rd_mode_ == 1 ? arma::dot(res, weights_) : res.max();
  }
  
  void optim(){
    arma::uword s_full_ = s_;
    new_val_ = comp_c_obj_fun();
    int i = 0;
    double diff = 1.0;
    double r = 0.01;
    while(diff > 0){
      i++;
      val_ = new_val_;
      if (trace_) Rcpp::Rcout << "\nIter " << i << ":Var: " << val_ << "\n";
      // make s swaps
      track_swaps_.fill(0);
      arma::mat unique_designs(k_,n_,fill::zeros);
      arma::vec track_swaps_val(n_,fill::zeros);
      track_swaps_val(0) = val_;
      arma::uword i_s_;
      arma::uword search_depth;
      bool r_loop = true;
      arma::uword min_val;
      while(r_loop){
        if (trace_) Rcpp::Rcout << "\rr: " << r;
        i_s_ = 0;
        search_depth=0;
        while(search_depth <= s_){
          // find optimum swap
          arma::mat swap_vals(k_,k_,fill::zeros);
#pragma omp parallel for
          for(arma::uword j=0; j<k_; j++){
            arma::vec val_in_vec = eval_swap(j);
            swap_vals.row(j) = val_in_vec.t();
          }
          bool idx_isin = true;
          arma::uword out_idx_sort = 0;
          //arma::uword out_idx_sort_col = 0;
          arma::uvec swap_vals_sort = sort_index(swap_vals);
          arma::uword target;
          arma::uword rm_target;
          double old_total;
          double des_total;
          while(idx_isin){
            target = floor((swap_vals_sort(out_idx_sort))/k_); //swap_idx(swap_vals_sort(out_idx_sort),out_idx_sort_col);
            rm_target = (swap_vals_sort(out_idx_sort)) - target*k_; //swap_vals_sort(out_idx_sort);
            unique_designs.col(i_s_+1) = unique_designs.col(i_s_);
            unique_designs(target,i_s_+1)++;
            unique_designs(rm_target,i_s_+1)--;
            old_total = sum(abs(unique_designs.col(i_s_)))/2;
            des_total = sum(abs(unique_designs.col(i_s_+1)))/2;
            //if (trace_) Rcpp::Rcout << "\ndes total " << des_total << " val: " << swap_vals(rm_target,target)  << target << rm_target;
            
            track_swaps_val(i_s_+1) = swap_vals(rm_target,target)+r*des_total;
            if(all(track_swaps_val.subvec(0,i_s_) > track_swaps_val(i_s_+1)) || 
               (any(track_swaps_val.subvec(0,i_s_) <= track_swaps_val(i_s_+1)) && des_total > old_total)){
              idx_isin = false;
            } else {
              out_idx_sort++;
            }
          }
          //if (trace_) Rcpp::Rcout << "\ni_s_: " << i_s_ << " depth: " << search_depth << " val: " << track_swaps_val(i_s_+1);
          
          //track_swaps_val(i_s_) = swap_vals(rm_target,target); 
          track_swaps_(i_s_,0) = rm_target;
          track_swaps_(i_s_,1) = target;
          make_swap(rm_target,target);
          Update_M_list();
          i_s_++;
          search_depth = des_total;
        }
        if(all(track_swaps_val.subvec(1,i_s_-1) >= val_) && r>1e-15){
          r *= 0.5;
        } else {
          min_val = index_min(track_swaps_val.subvec(1,i_s_));
          r_loop = false;
        }
      }
      
      if (trace_) Rcpp::Rcout << "obj: " << track_swaps_val.subvec(0,i_s_-1).t();
      if(any(track_swaps_val.subvec(1,i_s_-1) < val_)){
        if(min_val < (i_s_-1)){
          for(arma::uword is=i_s_;is>=min_val+2;is--){
            //if (trace_) Rcpp::Rcout << "\nBacktracking..." << is-1 << " of " << min_val;
            //arma::uvec rm_cond = find(idx_in_sub_==track_swaps_(i_s_,1));
            make_swap(track_swaps_(is-1,1),track_swaps_(is-1,0));
          }
        }
        Update_M_list();
        new_val_ = comp_c_obj_fun();
        diff = val_- new_val_ ;
        if (trace_)Rcpp::Rcout << "\nDiff " << diff;
        if(diff > 0){
          idx_in_ = idx_in_sub_;
          A_list_ = A_list_sub_;
          X_list_ = X_list_sub_;
          M_list_ = M_list_sub_;
          s_ = s_full_;
        }
      } else {
        diff = 0;
        //track_swaps_val.t().raw_print();
      }
    }
    Rcpp::Rcout << "\ncompleted, var: " << val_ << "v" << new_val_<< std::endl;
  }
  
private:
  
  void make_swap(arma::uword outobs, arma::uword inobs){
    curr_obs_(outobs)--;
    curr_obs_(inobs)++;
    // if(trace_){
    //   bool outbool = any(idx_in_sub_ == outobs);
    //   Rcpp::Rcout<< "\nOut: " << outbool;
    // }
    arma::uvec rm_cond = find(idx_in_sub_== outobs);
    arma::uword rm_idx = rm_cond(0);
    idx_in_rm_ = uvec_minus(idx_in_sub_,rm_idx);
    idx_in_sub_ = join_idx(idx_in_rm_,inobs);
    for (arma::uword idx = 0; idx < nlist_; ++idx) {
      arma::mat A1 = A_list_sub_.slice(idx);
      const arma::mat rm1A = remove_one_many_mat(A1, arma::uvec({rm_idx}));
      arma::rowvec z_j = Z_all_list_(idx,0).row(inobs);
      arma::mat z_d = Z_all_list_(idx,0).rows(idx_in_rm_);
      double sig_jj = arma::as_scalar(z_j * D_list_(idx,0) * z_j.t()); 
      sig_jj += W_all_diag_(inobs);
      arma::vec f = z_d * D_list_(idx,0) * z_j.t();
      A_list_sub_.slice(idx) = add_one_mat(rm1A, 
                        sig_jj,
                        f);
      
      arma::mat X = X_list_sub_(idx,0);
      X.rows(0,n_-2) = X_list_sub_(idx,0).rows(uvec_minus(all_idx_,rm_idx));
      X.row(n_-1) = X_all_list_(idx,0).row(inobs);
      X_list_sub_(idx,0) = X;
    }
  }
  
  void Update_M_list(bool useA = false, bool M_check = false) {
    for (arma::uword i = 0; i < nlist_; ++i) {
      arma::mat X1 = useA ? X_list_(i,0) : X_list_sub_(i,0);
      arma::mat A = useA ? A_list_.slice(i) : 
        A_list_sub_.slice(i);
      arma::mat M = X1.t() * A * X1;
      if(useA){
        M_list_(i,0) = M;
        M_list_sub_(i,0) = M;
      } else {
        M_list_sub_(i,0) = M;
      }
      if (!M_check && !M.is_sympd()) {
        arma::vec colsum = arma::trans(arma::sum(M, 0));
        arma::uvec colidx = arma::find(colsum == 0);
        Rcpp::Rcout << "ERROR MSG:\n"
                    << "Design " << i+1 << " has the following column(s) with zero colsum:\n"
                    << colidx.t() + 1 << std::endl;
        Rcpp::stop("M not positive semi-definite.");
      }
    }
  }
  
  void UpdateResult(arma::uword idx_to_rm = 0,
                    arma::uword idx_to_add = 0) {
    // if we are keeping a swap then this function is called... 
    idx_in_ = idx_in_sub_;
    A_list_ = A_list_sub_;
    X_list_ = X_list_sub_;
    M_list_ = M_list_sub_;
  }
  
  
  arma::vec eval_swap(arma::uword obs){
    arma::mat val_in_mat(k_,nlist_,arma::fill::zeros);
    bool obsisin = any(idx_in_sub_ == obs);
    if(obsisin){
      arma::uvec rm_cond = find(idx_in_sub_== obs);
      arma::uword rm_idx = rm_cond(0);
      for (arma::uword j = 0; j < nlist_; ++j) {
        arma::mat A1 = A_list_sub_.slice(j);
        const arma::mat rm1A = remove_one_many_mat(A1, arma::uvec({rm_idx}));
        arma::mat X(n_,X_list_sub_(j,0).n_cols);
        X.rows(0,n_-2) = X_list_sub_(j,0).rows(uvec_minus(all_idx_,rm_idx));
        for (arma::uword i = 0; i < k_; ++i) {
          if(obs != i && curr_obs_(i)<max_obs_(i)){
            X.row(n_-1) = X_all_list_(j,0).row(i);
            arma::rowvec z_j = Z_all_list_(j,0).row(i);
            arma::mat z_d = Z_all_list_(j,0).rows(uvec_minus(idx_in_sub_,rm_idx));
            double sig_jj = arma::as_scalar(z_j * D_list_(j,0) * z_j.t()); 
            sig_jj += W_all_diag_(i);
            arma::vec f = z_d * D_list_(j,0) * z_j.t();
            arma::mat sub1 = add_one_mat(rm1A, 
                                         sig_jj,
                                         f);
            val_in_mat(i,j) = c_obj_fun(
              X.t() * sub1 * X,
              C_list_(j,0));
            // val_in_mat(i,j) = c_d_deriv(arma::conv_to<arma::mat>::from(M_list_(j,0)), 
            //           X.t() * sub1 * X, 
            //           C_list_(j,0));
          } else {
            val_in_mat(i,j) = 10000;
          }
        }
      }
    } else {
      val_in_mat.fill(10000);
    }
    
    return rd_mode_==1 ? val_in_mat * weights_ : arma::vec(arma::max(val_in_mat, 1));
  }
  
};


//' Hill-Climbing algorithm to identify optimal GLMM design
//' 
//' Hill-Climbing algorithm to identify optimal GLMM design
//' @param N Integer specifying number of experimental conditions in the optimal design
//' @param idx_in Integer vector specifying the indexes of the experimental conditions to start from
//' @param C_list List of C vectors for the c-optimal function, see \link{glmmr}[DesignSpace]
//' @param X_list List of X matrices
//' @param sig_list List of inverse covariance matrices
//' @param weights Vector specifying the weights of each design
//' @param exp_cond Vector specifying the experimental condition index of each observation
//' @param nfix Vector listing the experimental condition indexes that are fixed in the design
//' @param any_fix Integer. 0 = no experimental conditions are fixed, 1 = some experimental conditions are fixed
//' @param rd_mode Integer. Robust objective function, 1=weighted average, 2=minimax
//' @param trace Logical indicating whether to provide detailed output
//' @return A vector of experimental condition indexes in the optimal design
// [[Rcpp::export]]
Rcpp::List GradRobustStep(arma::uvec idx_in, 
                          Rcpp::List C_list, 
                          Rcpp::List X_list, 
                          Rcpp::List Z_list, 
                          Rcpp::List D_list, 
                          arma::mat w_diag,
                          arma::uvec max_obs,
                          arma::vec weights,
                          arma::uvec nfix, 
                          arma::uword s=10,
                          double r = 0.01,
                          double b = 1,
                          arma::uword any_fix = 0,
                          arma::uword rd_mode = 1,
                          bool trace = true) {
  // need to map Rcpp list to the field here:
  arma::uword ndesign = weights.n_elem;
  arma::field<arma::vec> Cfield(ndesign);
  arma::field<arma::mat> Xfield(ndesign);
  arma::field<arma::mat> Zfield(ndesign);
  arma::field<arma::mat> Dfield(ndesign);
  for(arma::uword j=0; j<ndesign; j++){
    Cfield[j] = as<arma::vec>(C_list[j]);
    Xfield[j] = as<arma::mat>(X_list[j]);
    Zfield[j] = as<arma::mat>(Z_list[j]);
    Dfield[j] = as<arma::mat>(D_list[j]);
  }
  HillClimbing hc(idx_in, Cfield, Xfield, Zfield, Dfield,
                  w_diag,max_obs,
                  weights, any_fix, nfix,s,r,b, rd_mode, trace);
  hc.optim();
  return Rcpp::List::create(Named("idx_in") = hc.idx_in_,
                            Named("best_val_vec") = hc.best_val_vec_);
}
