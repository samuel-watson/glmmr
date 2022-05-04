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
  const arma::field<arma::mat> D_list_;
  const arma::field<arma::mat> X_all_list_; 
  const arma::field<arma::mat> Z_all_list_; 
  const arma::mat W_all_diag_;
  const arma::vec weights_; // weights for each design
  arma::uvec max_obs_;
  arma::uvec curr_obs_;
  const arma::uword nlist_; // number of designs
  arma::uword any_fix_;
  arma::uword n_; //size of the design to find
  arma::uword k_; //unique number of experimental conditions
  arma::uword nmax_;
  arma::uword p_;
  arma::uword q_;
  
public:
  arma::uvec idx_in_; 
  arma::uvec idx_in_sub_;
  arma::uvec idx_in_rm_;
  const arma::uvec exp_cond_; //vector linking row to experimental condition number
  arma::uword r_in_design_;
  arma::uword r_in_rm_;
  arma::uvec rows_in_design_;
  arma::uvec count_exp_cond_;
  arma::uvec count_exp_cond_rm_;
  double val_; // overall best value
  double new_val_; // new value
  double rm_val_;
  arma::uword fcalls_;
  arma::uword matops_;
  
  arma::cube A_list_; // inverse sigma matrices
  arma::cube rm1A_list_; // inverse sigma matrices with one removed - initialised to minus one but now needs to resize
  arma::cube M_list_;
  arma::cube M_list_sub_;
  const arma::uvec nfix_; //the indexes of the experimental conditions to keep
  const arma::uword rd_mode_; // robust designs mode: 1 == weighted, 2 == minimax.
  
  bool trace_;
  bool uncorr_;
  
public:
  HillClimbing(arma::uvec idx_in, 
               arma::uword n,
               arma::field<arma::vec> C_list, 
               arma::field<arma::mat> X_list, 
               arma::field<arma::mat> Z_list, 
               arma::field<arma::mat> D_list,
               arma::mat w_diag,
               arma::uvec max_obs,
               arma::vec weights,
               arma::uvec exp_cond,
               arma::uword any_fix,
               arma::uvec nfix, 
               arma::uword rd_mode = 0, 
               bool trace=false,
               bool uncorr=false) :
  C_list_(C_list), 
  D_list_(D_list),
  X_all_list_(X_list),
  Z_all_list_(Z_list),
  W_all_diag_(w_diag),
  weights_(weights), 
  max_obs_(max_obs),
  curr_obs_(max_obs.n_elem,fill::zeros),
  nlist_(weights_.n_elem),
  any_fix_(any_fix),
  n_(n), 
  k_(max_obs.n_elem),
  nmax_(2*ceil(X_all_list_(0,0).n_rows/k_)*n_),
  p_(X_all_list_(0,0).n_cols),
  q_(Z_all_list_(0,0).n_cols),
  idx_in_(idx_in),
  idx_in_sub_(idx_in),
  idx_in_rm_(idx_in),
  exp_cond_(exp_cond),
  r_in_design_(0),
  r_in_rm_(0),
  rows_in_design_(nmax_,arma::fill::zeros),
  count_exp_cond_(nmax_,arma::fill::zeros),
  count_exp_cond_rm_(nmax_,arma::fill::zeros),
  val_(0.0), 
  new_val_(0.0),
  rm_val_(0.0),
  fcalls_(0),
  matops_(0),
  A_list_(nmax_,nmax_,nlist_,fill::zeros), 
  rm1A_list_(nmax_,nmax_,nlist_,fill::zeros),
  M_list_(p_,p_,nlist_,fill::zeros),
  M_list_sub_(p_,p_,nlist_,fill::zeros),
  nfix_(nfix), 
  rd_mode_(rd_mode), 
  trace_(trace),
  uncorr_(uncorr){
    build_XZ();
  }
  
  arma::uvec join_idx(const arma::uvec &idx, arma::uword elem) {
    return arma::join_cols(idx, arma::uvec({elem}));
  }
  
  void build_XZ(){
    //Rcpp::Rcout << "\nStart XZ";
    idx_in_ = sort(idx_in_);
    for(arma::uword i=0; i<idx_in_.n_elem; i++){
      curr_obs_(idx_in_(i)-1)++;
    }
    idx_in_sub_ = idx_in_;
    arma::vec idx = arma::linspace(0, nmax_-1, nmax_);
    rows_in_design_ = arma::conv_to<arma::uvec>::from(idx);
    arma::vec vals(nlist_);
    for(arma::uword j=0; j<nlist_;j++){
      arma::mat X(nmax_,p_,fill::zeros);
      arma::mat Z(nmax_,q_,fill::zeros);
      arma::vec w_diag(nmax_);
      arma::uword rowcount = 0;
      for(arma::uword k=0; k<idx_in_.n_elem;k++){
        arma::uvec rowstoincl = find(exp_cond_ == idx_in_(k));
        
        //Rcpp::Rcout << "\nR incl: " << rowstoincl.t() << " X dim " << X_all_list_(j,0).n_rows << " " << X.n_rows<< " " << rowcount;
        for(arma::uword l=0;l<rowstoincl.n_elem;l++){
          X.row(rowcount) = X_all_list_(j,0).row(rowstoincl(l));
          Z.row(rowcount) = Z_all_list_(j,0).row(rowstoincl(l));
          w_diag(rowcount) = W_all_diag_(rowstoincl(l),j);
          count_exp_cond_(k)++;
          rowcount++;
        }
      }
      //Rcpp::Rcout << "\nFinish XZ";
      if(j==0)r_in_design_ = rowcount;
      arma::mat tmp = Z.head_rows(rowcount)*D_list_(j,0)*Z.head_rows(rowcount).t();
      tmp.diag() += w_diag.head(rowcount);
      
      //Rcpp::Rcout << "\nCount: " << rowcount;
      if(uncorr_){
        M_list_.slice(j) = X.head_rows(rowcount).t() * tmp.i() * X.head_rows(rowcount);
        M_list_sub_.slice(j) = M_list_.slice(j);
        vals(j) = c_obj_fun( M_list_.slice(j), C_list_(j,0));
      } else {
        A_list_.slice(j).submat(0,0,r_in_design_-1,r_in_design_-1) = tmp.i();
        vals(j) = c_obj_fun( X.head_rows(rowcount).t() * A_list_.slice(j).submat(0,0,r_in_design_-1,r_in_design_-1) * X.head_rows(rowcount), C_list_(j,0));
      }
      
    }
    new_val_ = rd_mode_ == 1 ? 1/arma::dot(vals, weights_) : 1/vals.max();
    //Rcpp::Rcout << "\nval: " << new_val_;
  }
  
  void local_search(){
    if (trace_) Rcpp::Rcout << "\nLocal search";
    int i = 0;
    double diff = 1.0;
    while(diff > 0){
      i++;
      val_ = new_val_;
      if (trace_) Rcpp::Rcout << "\nIter " << i << ": Prec: " << val_ << " Var: " << 1/val_;
      // evaluate the swaps
      arma::mat val_swap(k_,k_,fill::zeros);
      for(arma::uword j=1; j < k_+1; j++){
        if(any(idx_in_ == j)){
          arma::vec val_in_vec = eval(true,j);
          val_swap.row(j-1) = val_in_vec.t();
        }
      }
      
      double newval = val_swap.max();
      diff = newval - val_;
      if(diff > 0){
        arma::uword swap_sort = val_swap.index_max();
        arma::uword target = floor(swap_sort/k_); 
        arma::uword rm_target = (swap_sort) - target*k_;
        if(uncorr_){
          rm_obs_uncor(rm_target+1);
          new_val_ = add_obs_uncor(target+1,true,true);
        } else {
          rm_obs(rm_target+1);
          new_val_ = add_obs(target+1,true,true);
        }
      }
    }
  }
  
  void greedy_search(){
    // step 1: find optimal smallest design
    int i = 0;
    if (trace_) Rcpp::Rcout << "\nidx: " << idx_in_.t();
    if (trace_) Rcpp::Rcout << "\nGreedy search: " << n_;
    arma::uword idxcount = idx_in_.n_elem;
    while(idxcount < n_){
      i++;
      idxcount++;
      val_ = new_val_;
      if (trace_) Rcpp::Rcout << "\nIter " << i << " size: " << idxcount << " Prec: " << val_ << " Var: " << 1/val_;
      arma::vec val_swap = eval(false);
      arma::uword swap_sort = val_swap.index_max();
      if (trace_) Rcpp::Rcout << " adding " << swap_sort+1;
      if(uncorr_){
        new_val_ = add_obs_uncor(swap_sort+1,false,true);
      } else {
        new_val_ = add_obs(swap_sort+1,false,true); 
      }
    }
  }
  
private:
  arma::uvec get_rows(arma::uword idx){
    arma::uword start = (idx == 0) ? 0 : sum(count_exp_cond_.subvec(0,idx-1));
    arma::uword end = start + count_exp_cond_(idx) - 1;
    return rows_in_design_.subvec(start,end);
  }
  
  arma::uvec get_all_rows(arma::uvec idx){
    // get the rows corresponding to the given vector of experimental conditions
    arma::uvec rowidx(nmax_);
    arma::uword count = 0;
    for(arma::uword i = 0; i < idx.n_elem; i++){
      arma::uvec addidx = find(exp_cond_ == idx(i));
      for(arma::uword j = 0; j < addidx.n_elem; j++){
        rowidx(count) = addidx(j);
        count++;
      }
    }
    return rowidx.head(count);
  }
  
  arma::uvec idx_minus(arma::uword idx){
    arma::uword start = (idx == 0) ? 0 : sum(count_exp_cond_.subvec(0,idx-1));
    arma::uword end = start + count_exp_cond_(idx) - 1;
    arma::uvec idxrtn(r_in_design_-(end-start+1));
    idxrtn.head(start) = rows_in_design_.head(start);
    if(end < r_in_design_-1)idxrtn.tail(r_in_design_-end) = rows_in_design_.subvec(end+1,r_in_design_-1);
    return rows_in_design_.subvec(start,end);
  }
  
  void rm_obs(arma::uword outobs){
    arma::uvec rm_cond = find(idx_in_== outobs);
    //Rcpp::Rcout << "\nRm: " << rm_cond.t();
    arma::uvec rowstorm = get_rows(rm_cond(0));
    for (arma::uword idx = 0; idx < nlist_; ++idx) {
      matops_++;
      arma::mat A1 = A_list_.slice(idx).submat(0,0,r_in_design_-1,r_in_design_-1);
      const arma::mat rm1A = remove_one_many_mat(A1, rowstorm);
      if(idx==0)r_in_rm_ = rm1A.n_rows;
      rm1A_list_.slice(idx).submat(0,0,r_in_rm_-1,r_in_rm_-1) = rm1A;
    }
    idx_in_rm_ = uvec_minus(idx_in_,rm_cond(0));
    count_exp_cond_rm_.head(rm_cond(0)) = count_exp_cond_.head(rm_cond(0));
    if(rm_cond(0)>=idx_in_.n_elem - 1){
      count_exp_cond_rm_(rm_cond(0)) = count_exp_cond_(rm_cond(0)+1);
    } else {
      count_exp_cond_rm_.subvec(rm_cond(0),idx_in_.n_elem-2) = count_exp_cond_.subvec(rm_cond(0)+1,idx_in_.n_elem-1);
    }
  }
  
  void rm_obs_uncor(arma::uword outobs){
    arma::uvec rm_cond = find(idx_in_ == outobs);
    arma::uvec rowstorm = find(exp_cond_ == outobs);
    for(arma::uword j=0; j<nlist_;j++){
      arma::mat X(rowstorm.n_elem,p_,fill::zeros);
      arma::mat Z(rowstorm.n_elem,q_,fill::zeros);
      arma::vec w_diag(rowstorm.n_elem);
      
      for(arma::uword l=0;l<rowstorm.n_elem;l++){
        X.row(l) = X_all_list_(j,0).row(rowstorm(l));
        Z.row(l) = Z_all_list_(j,0).row(rowstorm(l));
        w_diag(l) = W_all_diag_(rowstorm(l),j);
      }
      
      arma::mat tmp = Z*D_list_(j,0)*Z.t();
      tmp.diag() += w_diag;
      M_list_sub_.slice(j) = M_list_.slice(j) - X.t() * tmp.i() * X;
    }
    idx_in_rm_ = uvec_minus(idx_in_,rm_cond(0));
    count_exp_cond_rm_.head(rm_cond(0)) = count_exp_cond_.head(rm_cond(0));
    if(rm_cond(0)>=idx_in_.n_elem - 1){
      count_exp_cond_rm_(rm_cond(0)) = count_exp_cond_(rm_cond(0)+1);
    } else {
      count_exp_cond_rm_.subvec(rm_cond(0),idx_in_.n_elem-2) = count_exp_cond_.subvec(rm_cond(0)+1,idx_in_.n_elem-1);
    }
  }
  
  double add_obs(arma::uword inobs,
                 bool userm = true,
                 bool keep = false){
    arma::uvec rowstoadd = find(exp_cond_ == inobs);
    arma::uvec idxvec = userm ? idx_in_rm_ : idx_in_;
    arma::uvec idxexist = get_all_rows(idxvec);
    arma::uword n_to_add = rowstoadd.n_elem;
    arma::uword n_already_in = idxexist.n_elem;
    arma::uvec idx_in_vec(n_already_in + n_to_add, fill::zeros);
    idx_in_vec(span(0,n_already_in - 1)) = idxexist;
    arma::mat X(n_already_in + n_to_add,p_,fill::zeros);
    arma::vec vals(nlist_);
    bool issympd = true;
    for (arma::uword idx = 0; idx < nlist_; ++idx) {
      arma::mat A = userm ? rm1A_list_.slice(idx).submat(0,0,r_in_rm_-1,r_in_rm_-1) : 
        A_list_.slice(idx).submat(0,0,r_in_design_-1,r_in_design_-1);
      X.rows(span(0,n_already_in-1)) = X_all_list_(idx,0).rows(idxexist);
      for(arma::uword j = 0; j < n_to_add; j++){
        arma::rowvec z_j = Z_all_list_(idx,0).row(rowstoadd(j));
        arma::mat z_d = Z_all_list_(idx,0).rows(idx_in_vec(span(0,n_already_in - 1)));
        double sig_jj = arma::as_scalar(z_j * D_list_(idx,0) * z_j.t()); 
        sig_jj += W_all_diag_(rowstoadd(j));
        arma::vec f = z_d * D_list_(idx,0) * z_j.t();
        A = add_one_mat(A, 
                        sig_jj,
                        f);
        idx_in_vec(n_already_in) = rowstoadd(j);
        X.row(n_already_in) = X_all_list_(idx,0).row(rowstoadd(j));
        n_already_in++;
      }
      //check if positive definite
      arma::mat M = X.t() * A * X;
      issympd = M.is_sympd();
      if(issympd){
        if(keep){
          if(idx==0)r_in_design_ = A.n_rows;
          A_list_.slice(idx).submat(0,0,r_in_design_-1,r_in_design_-1) = A;
        }
        vals(idx) = c_obj_fun( M, C_list_(idx,0));
      } else {
        vals.fill(0);
        break;
      }
    }
    if(keep && issympd){
      if(userm){
        idx_in_ = join_idx(idx_in_rm_,inobs);
        curr_obs_(inobs-1)++;
        count_exp_cond_.subvec(0,idx_in_.n_elem-2) = count_exp_cond_rm_.subvec(0,idx_in_.n_elem-2);
        count_exp_cond_(idx_in_.n_elem-1) = n_to_add;
      } else {
        idx_in_ = join_idx(idx_in_,inobs);
        curr_obs_(inobs-1)++;
        count_exp_cond_(idx_in_.n_elem-1) = n_to_add;
      }
    }
    double rtn = rd_mode_ == 1 ? arma::dot(vals, weights_) : vals.max();
    if(rtn > 0){
      return 1/rtn;
    } else {
      return 0;
    }
  }
  
  double add_obs_uncor(arma::uword inobs,
                       bool userm = true,
                       bool keep = false){
    arma::vec vals(nlist_);
    arma::uvec rowstoadd = find(exp_cond_ == inobs);
    bool issympd = true;
    for(arma::uword j=0; j<nlist_;j++){
      arma::mat X(rowstoadd.n_elem,p_,fill::zeros);
      arma::mat Z(rowstoadd.n_elem,q_,fill::zeros);
      arma::vec w_diag(rowstoadd.n_elem);
      
      for(arma::uword l=0;l<rowstoadd.n_elem;l++){
        X.row(l) = X_all_list_(j,0).row(rowstoadd(l));
        Z.row(l) = Z_all_list_(j,0).row(rowstoadd(l));
        w_diag(l) = W_all_diag_(rowstoadd(l),j);
      }
      
      arma::mat tmp = Z*D_list_(j,0)*Z.t();
      tmp.diag() += w_diag;
      arma::mat M = userm ? M_list_sub_.slice(j) : M_list_.slice(j);
      
      M += X.t() * tmp.i() * X;
      //Rcpp::Rcout << "\nid:" << inobs << " M: " << M;
      issympd = M.is_sympd();
      if(issympd){
        if(keep){
          M_list_.slice(j) = M;
        }
        vals(j) = c_obj_fun(M, C_list_(j,0));
      } else {
        vals.fill(0);
        break;
      }
    }
    if(keep && issympd){
      if(userm){
        idx_in_ = join_idx(idx_in_rm_,inobs);
        curr_obs_(inobs-1)++;
        count_exp_cond_.subvec(0,idx_in_.n_elem-2) = count_exp_cond_rm_.subvec(0,idx_in_.n_elem-2);
        count_exp_cond_(idx_in_.n_elem-1) = rowstoadd.n_elem;
      } else {
        idx_in_ = join_idx(idx_in_,inobs);
        curr_obs_(inobs-1)++;
        count_exp_cond_(idx_in_.n_elem-1) = rowstoadd.n_elem;
      }
    }
    double rtn = rd_mode_ == 1 ? arma::dot(vals, weights_) : vals.max();
    if(rtn > 0){
      return 1/rtn;
    } else {
      return 0;
    }
  }
  
  arma::vec eval(bool userm = true, arma::uword obs = 0){
    arma::vec val_in_mat(k_,arma::fill::zeros);
    if(userm){
      bool obsisin = any(idx_in_ == obs);
      if(obsisin){
        if(uncorr_){
          rm_obs_uncor(obs);
        } else {
          rm_obs(obs);
        }
#pragma omp parallel for
        for (arma::uword i = 1; i < k_+1; ++i) {
          if(obs != i && curr_obs_(i-1)<max_obs_(i-1)){
            if(uncorr_){
              val_in_mat(i-1) = add_obs_uncor(i,true,false);
              
            } else {
              val_in_mat(i-1) = add_obs(i,true,false);
            }
          } 
        }
        matops_ += k_*nlist_;
        fcalls_ += k_*nlist_;
      } 
    } else {
#pragma omp parallel for
      for (arma::uword i = 1; i < k_+1; ++i) {
        if(curr_obs_(i-1)<max_obs_(i-1)){
          if(uncorr_){
            val_in_mat(i-1) = add_obs_uncor(i,false,false);
          } else {
            val_in_mat(i-1) = add_obs(i,false,false);
          }
        }
      }
      matops_ += k_*nlist_;
      fcalls_ += k_*nlist_;
    }
    
    return val_in_mat;
  }
  
};


//' Hill-Climbing algorithm to identify optimal GLMM design
//' 
//' Hill-Climbing algorithm to identify optimal GLMM design
//' @param N Integer specifying number of experimental conditions in the optimal design
//' @param idx_in Integer vector specifying the indexes of the experimental conditions to start from
//' @param n Integer specifying the size of the design to find. For local search, this should be equal to the size of idx_in
//' @param C_list List of C vectors for the c-optimal function, see \link{glmmr}[DesignSpace]
//' @param X_list List of X matrices
//' @param sig_list List of inverse covariance matrices
//' @param weights Vector specifying the weights of each design
//' @param exp_cond Vector specifying the experimental condition index of each observation
//' @param nfix Vector listing the experimental condition indexes that are fixed in the design
//' @param any_fix Integer. 0 = no experimental conditions are fixed, 1 = some experimental conditions are fixed
//' @param type Integer. 0 = local search algorith. 1 = greedy search algorithm.
//' @param rd_mode Integer. Robust objective function, 1=weighted average, 2=minimax
//' @param trace Logical indicating whether to provide detailed output
//' @return A vector of experimental condition indexes in the optimal design
// [[Rcpp::export]]
Rcpp::List GradRobustStep(arma::uvec idx_in, 
                          arma::uword n,
                          Rcpp::List C_list, 
                          Rcpp::List X_list, 
                          Rcpp::List Z_list, 
                          Rcpp::List D_list, 
                          arma::mat w_diag,
                          arma::uvec max_obs,
                          arma::vec weights,
                          arma::uvec exp_cond,
                          arma::uvec nfix, 
                          arma::uword any_fix = 0,
                          arma::uword type = 0,
                          arma::uword rd_mode = 1,
                          bool trace = true,
                          bool uncorr = false) {
  // need to map Rcpp list to the field here:
  arma::uword ndesign = weights.n_elem;
  arma::field<arma::vec> Cfield(ndesign);
  arma::field<arma::mat> Xfield(ndesign);
  arma::field<arma::mat> Zfield(ndesign);
  arma::field<arma::mat> Dfield(ndesign);
  for(arma::uword j=0; j<ndesign; j++){
    Cfield(j,0) = as<arma::vec>(C_list[j]);
    Xfield(j,0) = as<arma::mat>(X_list[j]);
    Zfield(j,0) = as<arma::mat>(Z_list[j]);
    Dfield(j,0) = as<arma::mat>(D_list[j]);
  }
  HillClimbing hc(idx_in, n, Cfield, Xfield, Zfield, Dfield,
                  w_diag,max_obs,
                  weights, exp_cond, any_fix, nfix,rd_mode, trace, uncorr);
  if(type==0)hc.local_search();
  if(type==1){
    hc.local_search();
    hc.greedy_search();
    hc.local_search();
  }
  if(type==2){
    hc.local_search();
    hc.greedy_search();
  }
  if(type==3){
    hc.greedy_search();
    hc.local_search();
  }
  return Rcpp::List::create(Named("idx_in") = hc.idx_in_,
                            Named("best_val_vec") = hc.val_,
                            Named("func_calls") = hc.fcalls_,
                            Named("mat_ops") = hc.matops_);
}
