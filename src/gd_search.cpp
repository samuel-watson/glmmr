#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;

#ifdef _OPENMP
#include <omp.h>
#endif

// [[Rcpp::depends(RcppArmadillo)]]

// [[Rcpp::export]]
double obj_fun(const arma::mat &A, const arma::vec &U2) {
  // this is the directional derivative
  return arma::as_scalar(U2.t() * A * U2);
}

// [[Rcpp::export]]
double c_obj_fun(arma::mat M, arma::vec C) {
  // this is the objective function c-optimal
  arma::mat M_inv = arma::inv_sympd(M);
  return arma::as_scalar(C.t() * M_inv * C);
}

// [[Rcpp::export]]
arma::mat gen_m(const arma::mat &X, const arma::mat &A) {
  //generate information matrix
  return X.t() * A * X;
}

arma::uvec match_uvec(arma::uvec x, arma::uword val){
  arma::uvec match_vec(x.n_elem, fill::value(-1));
  arma::uword count = 0;
  for(arma::uword j=0; j<x.n_elem; j++){
    if(x(j)==val){
      match_vec(count) = j;
      count += 1;
    }
  }
  return match_vec.rows(0,count-1);
}

// [[Rcpp::export]]
arma::mat remove_one_many_mat(const arma::mat &A, 
                              const arma::uvec &i) {
  arma::vec idx = arma::linspace(0, A.n_rows - 1, A.n_rows);
  arma::uvec uidx = arma::conv_to<arma::uvec>::from(idx);
  arma::uvec isort = arma::sort(i,"descend"); // sort descending to avoid row conflicts
  arma::mat A2 = A;
  
  for(arma::uword j=0; j<i.n_elem; j++){
    arma::vec idx_new = arma::linspace(0, uidx.n_elem - 1, uidx.n_elem);
    arma::uvec uidx_new = arma::conv_to<arma::uvec>::from(idx_new);
    uidx.shed_row(isort(j));
    uidx_new.shed_row(isort(j));
    double d = A2(isort(j), isort(j));
    arma::vec b = A2.submat(uidx_new, arma::uvec({isort(j)}));
    A2 = A2.submat(uidx_new, uidx_new) - b * b.t() / d;
  }
  
  return A2;
}

// [[Rcpp::export]]
double remove_one_many(const arma::mat &A, 
                       const arma::uvec &i,
                       const arma::vec &u) {
  arma::vec idx = arma::linspace(0, A.n_rows - 1, A.n_rows);
  arma::uvec uidx = arma::conv_to<arma::uvec>::from(idx);
  arma::uvec isort = arma::sort(i,"descend"); // sort descending to avoid row conflicts
  arma::mat A2 = A;
  for(arma::uword j=0; j<i.n_elem; j++){
    arma::vec idx_new = arma::linspace(0, uidx.n_elem - 1, uidx.n_elem);
    arma::uvec uidx_new = arma::conv_to<arma::uvec>::from(idx_new);
    uidx.shed_row(isort(j));
    uidx_new.shed_row(isort(j));
    double d = A2(isort(j), isort(j));
    arma::vec b = A2.submat(uidx_new, arma::uvec({isort(j)}));
    A2 = A2.submat(uidx_new, uidx_new) - b * b.t() / d;
  }
  return obj_fun(A2, u(uidx));
}

// [[Rcpp::export]]
double add_one(const arma::mat &A, 
               double sigma_jj, 
               const arma::vec &f,
               const arma::vec &u) {
  arma::mat A2(A.n_rows + 1, A.n_cols + 1, arma::fill::zeros);
  A2.submat(0, 0, A2.n_rows - 2, A2.n_cols - 2) = A;
  A2(A2.n_rows - 1, A2.n_cols - 1) = 1 / sigma_jj;
  
  // step 3: compute K2_inv
  arma::vec u1 = arma::join_cols(f, arma::vec({0}));
  arma::vec v1(u1.n_elem, arma::fill::zeros);
  v1(v1.n_elem - 1) = 1;
  A2 = A2 -
    ((A2 * u1) * (v1.t() * A2)) / (1 + arma::as_scalar((v1.t() * A2) * u1));
  
  // step 4: compute K3_inv
  A2 = A2 -
    ((A2 * v1) * (u1.t() * A2)) / (1 + arma::as_scalar((u1.t() * A2) * v1));
  
  return obj_fun(A2, u);
}

// [[Rcpp::export]]
arma::mat add_one_mat(const arma::mat &A, 
                      double sigma_jj, 
                      const arma::vec &f) {
  arma::mat A2(A.n_rows + 1, A.n_cols + 1, arma::fill::zeros);
  // step 1: compute A*
  A2.submat(0, 0, A2.n_rows - 2, A2.n_cols - 2) = A;
  for (arma::uword j = 0; j < A2.n_rows - 1; j++) {
    A2(j, A2.n_cols - 1) = 0;
    A2(A2.n_rows - 1, j) = 0;
  }
  A2(A2.n_rows - 1, A2.n_cols - 1) = 1 / sigma_jj;
  
  // step 3: compute K2_inv
  arma::vec u1 = arma::join_cols(f, arma::vec({0}));
  arma::vec v1(u1.n_elem, arma::fill::zeros);
  v1(v1.n_elem - 1) = 1;
  A2 = A2 -
    ((A2 * u1) * (v1.t() * A2)) / (1 + arma::as_scalar((v1.t() * A2) * u1));
  
  // step 4: compute K3_inv
  A2 = A2 -
    ((A2 * v1) * (u1.t() * A2)) / (1 + arma::as_scalar((u1.t() * A2) * v1));
  
  return A2;
}

arma::uvec std_setdiff(arma::uvec &x, arma::uvec &y) {
  std::vector<int> a = arma::conv_to<std::vector<int> >::from(arma::sort(x));
  std::vector<int> b = arma::conv_to<std::vector<int> >::from(arma::sort(y));
  std::vector<int> out;
  
  std::set_difference(a.begin(), a.end(), b.begin(), b.end(),
                      std::inserter(out, out.end()));
  
  return arma::conv_to<arma::uvec>::from(out);
}


// [[Rcpp::export]]
arma::uvec uvec_minus(const arma::uvec &v, arma::uword rm_idx) {
  arma::uword n = v.size();
  if (rm_idx == 0) return v.tail(n-1);
  if (rm_idx == n-1) return v.head(n-1);
  arma::uvec res(v.size()-1);
  res.head(rm_idx) = v.head(rm_idx);
  res.tail(n-1-rm_idx) = v.tail(n-1-rm_idx);
  return res;
}

class HillClimbing {
private:
  const arma::field<arma::vec> C_list_; // C vector
  const arma::field<arma::mat> X_list_; // cbind of X matrices
  const arma::cube sig_list_; //cbind of covariance matrices
  const arma::vec weights_; // weights for each design
  const arma::uword nlist_; // number of designs
  arma::uword any_fix_;
  arma::uword n_cond_; // number of experimental conditions
  const arma::uword n_; //size of the design to find
  const arma::uword N_; //size of the design space
  
  
public:
  
  //these are the indexes of the epxerimental conditions in and out of the design
  arma::uvec idx_in_; 
  arma::uvec idx_out_; 
  const arma::uvec exp_cond_; //vector linking row to experimental condition number 
  // these are the experimental conditions and row numbers relating to the design space of those currently in 
  // the design
  arma::uvec exp_cond_in_design_; 
  arma::uvec rows_in_design_; 
  arma::uword count_in_design_;
  // these are the experimental conditions and row numbers relating to the design space of those currently in 
  // rm1A 
  arma::uvec exp_cond_in_rm1_; // experimental conditions of the rows in the design
  arma::uvec rows_in_rm1_; //rows in current design relating to Sigma matrices
  arma::uword count_in_rm1_;
  // these are the experimental conditions and row numbers relating to the design space of those currently in 
  // a sub matrices
  arma::uvec exp_cond_in_asub_; // experimental conditions of the rows in the design
  arma::uvec rows_in_asub_; //rows in current design relating to Sigma matrices
  arma::uword count_in_asub_;
  
  arma::vec new_val_vec_; //value of objective function for each separate design
  arma::vec best_val_vec_; //value of objective function best
  double val_; // overall best value
  double new_val_; // new value
  arma::uvec ignore_idx_; // which to ignore if n_fix
  
  arma::cube A_list_; // cbind of inverse sigma matrices
  arma::cube A_list_sub_; // holder of 
  arma::cube rm1A_list_; // cbind of inverse sigma matrices with one removed - initialised to minus one but now needs to resize
  
  arma::field<arma::mat> M_list_; // field of information matrices
  arma::mat u_list_; // u = c^T M^-1 X^T rbind
  
  arma::uword A_nrows_; 
  arma::uword rm1A_nrows_;
  
  const arma::uvec nfix_; //CHANGED - this is now the indexes of the experimental conditions to keep
  
  //arma::uword n_nfix_; // number of fixed observations !! ADD IN 
  const arma::uword rd_mode_; // robust designs mode: 1 == weighted, 2 == minimax.
  
  bool trace_;
  
public:
  HillClimbing(arma::uvec idx_in, 
               arma::field<arma::vec> C_list, 
               arma::field<arma::mat> X_list, 
               arma::cube sig_list, 
               arma::vec weights,
               arma::uvec exp_cond,
               arma::uword any_fix,
               arma::uvec nfix, 
               arma::uword rd_mode = 0, 
               bool trace=false) :
  C_list_(C_list), 
  X_list_(X_list),
  sig_list_(sig_list),
  weights_(weights), 
  nlist_(weights_.n_elem),
  any_fix_(any_fix),
  n_cond_(),
  n_(idx_in.n_elem), 
  N_(sig_list_.n_rows),
  idx_in_(idx_in), 
  idx_out_(),
  exp_cond_(exp_cond),
  exp_cond_in_design_(N_),
  rows_in_design_(N_),
  count_in_design_(0),
  exp_cond_in_rm1_(N_),
  rows_in_rm1_(N_),
  count_in_rm1_(0),
  exp_cond_in_asub_(N_),
  rows_in_asub_(N_),
  count_in_asub_(0),
  new_val_vec_(nlist_, arma::fill::zeros), 
  best_val_vec_(nlist_, arma::fill::zeros), 
  val_(0.0), 
  new_val_(0.0),
  A_list_(N_,N_,nlist_,fill::zeros), 
  A_list_sub_(N_,N_,nlist_,fill::zeros),
  rm1A_list_(N_,N_,nlist_,fill::zeros),
  M_list_(nlist_,1),
  u_list_(N_, nlist_, arma::fill::zeros),
  A_nrows_(n_), 
  rm1A_nrows_(), 
  nfix_(nfix), 
  rd_mode_(rd_mode), 
  trace_(trace) {
    init_exp_cond();
    init_idx_out();
    Update_A_list();
    Update_M_list(0,true);
    Update_u_list();
    cout << "rd_mode = " << rd_mode_ << endl;
  }
  
  
  arma::uvec join_idx(const arma::uvec &idx, arma::uword elem) {
    return arma::join_cols(idx, arma::uvec({elem}));
  }
  
  // get the row numbers from exp_cond for a vector idx 
  // really a vectorised version of match - it matches all the idx in idx_list
  arma::uvec rows_to_keep(arma::uvec idx_list, 
                          arma::uvec idx){
    arma::uvec newrowstokeep(idx_list.n_elem,fill::zeros);
    arma::uword count_idx = 0;
    for(arma::uword j=0; j<idx.n_elem; j++){
      arma::uvec i_idx = match_uvec(idx_list,idx(j));
      newrowstokeep.rows(count_idx, count_idx+i_idx.n_elem-1) = i_idx;
      count_idx += i_idx.n_elem;
    }
    arma::vec newrowsvec = arma::conv_to<arma::vec>::from(newrowstokeep.subvec(0,count_idx-1));
    arma::uvec sortrows = arma::conv_to<arma::uvec>::from(arma::sort(newrowsvec));
    return sortrows;
  }
  
  void init_idx_out() {
    // calculate unique experimental conditions
    arma::uvec unique_exp_cond = unique(exp_cond_);
    n_cond_ = unique_exp_cond.n_elem;
    
    // generate the complete index
    arma::vec idx = arma::linspace(1, n_cond_, n_cond_);
    arma::uvec uidx = arma::conv_to<arma::uvec>::from(idx);
    // get the experimental conditions not included in idx_in_
    idx_out_ =  std_setdiff(uidx, idx_in_);
  }
  
  //initialise the list of indices currently in the design based on idx_in_
  void init_exp_cond(){
    arma::uvec inlist(N_,fill::zeros);
    arma::uvec rowlist(N_,fill::zeros);
    arma::uword count=0;
    
    for(arma::uword j=0; j<N_; j++){
      if(any(idx_in_ == exp_cond_(j))){
        inlist(count) = exp_cond_(j);
        rowlist(count) = j;
        count++;
      }
    }
    
    
    exp_cond_in_design_ = inlist;
    rows_in_design_ = rowlist;
    count_in_design_ = count;
    
    // also initialise the indexes to ignore
    if(any_fix_ > 0){
      arma::uvec ignore(n_,fill::zeros);
      arma::uword count_ignore = 0;
      for(arma::uword j=0; j<n_; j++){
        bool any_to_ignore = any(nfix_ == idx_in_(j));
        if(any_to_ignore){
          ignore(count_ignore) = j;
          ++count_ignore;
        }
      }
      ignore_idx_ = ignore.rows(0,count_ignore-1);
    }
  }
  
  //remove an experimental design from the list for rm1
  void rm_exp_cond(arma::uword idx){
    arma::uvec newexpcond = find(exp_cond_in_design_ != idx && exp_cond_in_design_ != 0);
    exp_cond_in_rm1_.fill(0);
    rows_in_rm1_.fill(0);
    count_in_rm1_ = newexpcond.n_elem;
    exp_cond_in_rm1_(span(0,count_in_rm1_ - 1)) = exp_cond_in_design_.elem(newexpcond);
    rows_in_rm1_(span(0,count_in_rm1_ - 1)) = rows_in_design_.elem(newexpcond);
  }
  
  //add an experimental condition at the bottom of the list
  // from rm1
  void add_exp_cond(arma::uword idx){
    arma::uvec newexpcondidx = find(exp_cond_ == idx);
    arma::uword n_add = newexpcondidx.n_elem;
    count_in_asub_ = count_in_rm1_ + n_add;
    exp_cond_in_asub_ = exp_cond_in_rm1_;
    rows_in_asub_ = rows_in_rm1_;
    rows_in_asub_(span(count_in_rm1_,count_in_asub_ - 1)) = newexpcondidx;
    for(arma::uword j=count_in_rm1_; j< count_in_asub_; j++){
      exp_cond_in_asub_(j) = idx;
    }
  }
  
  //index of the condition to remove next
  arma::uword idx_to_rm(){
    arma::uword idx = 0;
    if(any_fix_ > 0){
      for(arma::uword j=0; j< n_; j++){
        if(any(ignore_idx_ == j)){
          idx++;
        } else {
          break;
        }
      }
    } 
    return idx;
  }
  
  double comp_c_obj_fun() {
    for (arma::uword j = 0; j < nlist_; ++j) {
      new_val_vec_(j) = c_obj_fun(arma::conv_to<arma::mat>::from(M_list_(j,0)), 
                   C_list_(j,0));
    }
    return rd_mode_ == 1 ? arma::dot(new_val_vec_, weights_) : new_val_vec_.max();
  }
  
  double check_all_neighbours(double &diff) {
    // if reordering doesn't find a better solution then check all the neighbours
    // check optimality and see if there are any neighbours that improve the solution
    if (trace_) Rcpp::Rcout << "Checking optimality..." << std::endl;
    // looping through all elements of idx_in_
    for (arma::uword obs = 0; obs < n_; ++obs)  {
      if (trace_) Rcpp::Rcout << "\rChecking neighbour block: " << obs+1 << " of " << n_;
      bool check_skip = false;
      if(any_fix_ > 0){
        check_skip = any(ignore_idx_ == obs); // skip this block if it is in the fixed section
      }
      
      if(!check_skip){
        Update_rm1A_list(obs); // removes the current obs value 
        
        //bool idx_in_updated = false;
        
        // now loop through the idx_out_ elements to see if there are any winners
        // break the loop if we reach one with positive gradient
        for (arma::uword obs_j = 0; obs_j < idx_out_.n_elem; ++obs_j) {
          arma::vec val_in_mat(nlist_, arma::fill::zeros);
          // calculate the value of the observation
          for (arma::uword j = 0; j < nlist_; ++j) {
            val_in_mat(j) = add_one_helper_many(rm1A_list_.subcube(0,0,j,count_in_rm1_-1,count_in_rm1_-1,j), 
                       j, 
                       rows_in_rm1_(span(0,count_in_rm1_-1)), 
                       match_uvec(exp_cond_,idx_out_(obs_j)));
          }
          
          double val_in = rd_mode_==1 ? arma::sum(val_in_mat * weights_) : val_in_mat.max();
          // if the derivative is negative then go to next step, otherwise break this loop
          
          if (val_ - val_in < 0) {
            //add the observations to rm1 matrices
            Update_A_list_sub(obs_j);
            //generate new M matrices
            bool flag = Update_M_list(1);
            if (!flag) continue; // does this skip to the next part of the loop if not sym PD?
            // calculate new c-objective function
            new_val_ = comp_c_obj_fun();
            // if the variance is lower then keep the swap 
            if (new_val_ - val_ < 0) {
              diff = new_val_ - val_;
              if (trace_) Rcpp::Rcout << "\nImprovement found: " << new_val_ << std::endl;
              UpdateResult(obs, obs_j);
              //idx_in_updated = true;
              break;
            } // if (new_val - val < 0)
          } else break; // if (val - val_in < 0)
        } // for loop obs_j
        
        if (diff < 0) break;
      }
    } // for loop obs
    return diff;
  }
  
  // make a swap and compare values with or without reordering
  bool easy_swap(int &i, double &diff, bool reorder = false) {
    if (reorder) reorder_obs();
    arma::uword rm_idx = idx_to_rm(); //nfix - if nfix has length, then check which indexes to ignore
    Update_rm1A_list(rm_idx);  
    Update_A_list_sub(0);
    // if(trace_)Rcpp::Rcout << "DESIGN: " << rows_in_design_(span(0,count_in_design_-1)) << "\n ASUB: " <<
    //   rows_in_asub_(span(0,count_in_asub_-1));
    Update_M_list(2);
    new_val_ = comp_c_obj_fun();
    diff = new_val_ - val_;
    
    // we are now looking for the smallest value rather than largest so diff < 0
    std::string str = reorder ? "Iter (reorder)" : "Iter ";
    if (trace_) Rcpp::Rcout << str << i << ": " << val_ << std::endl;
    
    bool flag = diff < 0;
    if (flag) UpdateResult(rm_idx,0); // update the indexes of idx_in_ and idx_out_
    // if we don't get an improvement then we need to reset the obs in design
    
    return flag;
  }
  
  //MAIN ALGORITHM FUNCTION
  void grad_robust_step() {
    new_val_ = comp_c_obj_fun(); // value of objective function
    reorder_obs(); // order the experimental conditions in terms of their derivatives
    if(trace_)Rcpp::Rcout << "done reorder\n";
    double diff = -1.0;
    int i = 0;
    while (diff < 0) {
      ++i;
      if(trace_){
        arma::uword n_rep = n_ < 10 ? n_ : 10;
        idx_in_.head(n_rep).t().print("Exp. cond.: ");
      }
      val_ = new_val_; 
      best_val_vec_ = new_val_vec_;
      bool diff_l0 = easy_swap(i, diff); //returns true if an easy swap was made
      if (!diff_l0) {
        // if no easy swaps can be made, reorder the list
        diff_l0 = easy_swap(i, diff, true);
        // then check all the neighbours
        if (!diff_l0) check_all_neighbours(diff);
      }
    } // while loop
  }
  
  //OLD ALGORITHIM
  // UNTESTED AND DOESN'T REALLY WORK!
  void grad_robust_alg1() {
    new_val_ = comp_c_obj_fun();
    
    double diff = -1.0;
    int i = 0;
    // we now need diff to be negative
    while (diff < 0) {
      ++i;
      if(trace_){
        arma::uword n_rep = n_ < 10 ? n_ : 10;
        idx_in_.head(n_rep).t().print("Exp. cond.: ");
      }
      val_ = new_val_;
      reorder_obs();
      arma::uword rm_idx = idx_to_rm();
      Update_rm1A_list(rm_idx);  
      Update_A_list_sub(0);
      Update_M_list(2);
      new_val_ = comp_c_obj_fun();
      diff = new_val_ - val_;
      if (diff < 0) UpdateResult(rm_idx,0); 
      if (trace_) Rcpp::Rcout << "\rIter " << i << ": " << diff << " Var: " << new_val_ << std::endl;
    }
  }
  
private:
  
  // this just initialises the A _list
  void Update_A_list() {
    arma::uvec newrows = rows_in_design_(span(0,count_in_design_-1));
    for (arma::uword i = 0; i < nlist_; ++i) {
      const arma::mat tmp = sig_list_.slice(i).submat(newrows,newrows);
      A_list_.subcube(0,0,i,count_in_design_-1,count_in_design_-1,i) = tmp.i();
    }
    A_nrows_ = count_in_design_;
  }
  
  
  // adds observations to rm1a and puts in A_sub
  // ii is the index of idx_out_ to add
  void Update_A_list_sub(arma::uword ii) {
    //indexes to add using 
    add_exp_cond(idx_out_(ii));
    
    for(arma::uword i = 0; i < nlist_; ++i) {
      A_list_sub_.slice(i).fill(0);
      arma::mat rm1new = rm1A_list_.subcube(0,0,i,count_in_rm1_-1,count_in_rm1_-1,i);
      
      A_list_sub_.subcube(0,0,i,count_in_asub_-1,count_in_asub_-1,i) = add_one_helper_many_mat(rm1new, i,
                          rows_in_rm1_(span(0,count_in_rm1_-1)),
                          match_uvec(exp_cond_,idx_out_(ii)));
    }
    
  }
  
  // remove observations from A
  void Update_rm1A_list(arma::uword obs = 0) {
    rm_exp_cond(idx_in_(obs));
    arma::uvec i_idx = match_uvec(exp_cond_in_design_(span(0,count_in_design_-1)),idx_in_(obs));
    for (arma::uword idx = 0; idx < nlist_; ++idx) {
      rm1A_list_.slice(idx).fill(0);
      arma::mat A1 = A_list_.subcube(0,0,idx,count_in_design_-1,count_in_design_-1,idx);
      const arma::mat rm1A = remove_one_many_mat(A1, i_idx);
      rm1A_list_.subcube(0,0,idx,count_in_rm1_-1,count_in_rm1_-1,idx) = rm1A;
    }
    rm1A_nrows_ = rm1A_list_.n_rows; // update rows
  }
  
  
  // check_M == 0 ---- No positive definite check for M is needed
  //         == 1 ---- Do the check, return false if find not positive definite
  //         == 2 ---- Do the check, stop function if find not positive definite
  
  bool Update_M_list(int check_M = 0,
                     bool useA = false) {
    
    for (arma::uword i = 0; i < nlist_; ++i) {
      arma::mat X = X_list_(i,0);
      X = useA ? X.rows(rows_in_design_(span(0,count_in_design_-1))) :
        X.rows(rows_in_asub_(span(0,count_in_asub_-1))) ;
      
      arma::mat A = useA ? A_list_.subcube(0,0,i,count_in_design_-1,count_in_design_-1,i) : 
        A_list_sub_.subcube(0,0,i,count_in_asub_-1,count_in_asub_-1,i);
      arma::mat M = X.t() * A * X;
      M_list_(i,0) = M;
      if (check_M && !M.is_sympd()) {
        if (check_M == 1) return false;
        else {
          arma::vec colsum = arma::trans(arma::sum(M, 0));
          arma::uvec colidx = arma::find(colsum == 0);
          Rcpp::Rcout << "ERROR MSG:\n"
                      << "Design " << i+1 << " has the following column(s) with zero colsum:\n"
                      << colidx.t() + 1 << std::endl;
          Rcpp::stop("M not positive semi-definite.");
        }
      }
    }
    //M_list_ = newMlist;
    return true;
  }
  
  
  void Update_u_list() {
    for (arma::uword i = 0; i < nlist_; ++i) {
      arma::mat X = X_list_(i,0);
      arma::vec C = C_list_(i,0);
      arma::mat M_inv = arma::inv_sympd(arma::conv_to<arma::mat>::from(M_list_(i,0)));
      u_list_.col(i) = X * (M_inv * C);
    }
  }
  
  // if we are keeping a swap then this function is called...  
  void UpdateResult(arma::uword idx_to_rm = 0,
                    arma::uword idx_to_add = 0) {
    
    const arma::uword swap_idx = idx_in_(idx_to_rm); // exp cond being removed
    arma::uvec idx_in_shed = uvec_minus(idx_in_, idx_to_rm); 
    arma::uvec idx_out_shed = uvec_minus(idx_out_, idx_to_add);
    idx_in_ = join_idx(idx_in_shed, idx_out_(idx_to_add)); // combine the old idx_in_
    idx_out_ = join_idx(idx_out_shed, swap_idx);
    exp_cond_in_design_ = exp_cond_in_asub_;
    rows_in_design_ = rows_in_asub_;
    count_in_design_ = count_in_asub_;
    A_list_ = A_list_sub_;
    Update_u_list();
  }
  
  // evaluate remove_one
  
  arma::vec comp_val_out_vec() {
    arma::mat val_out_mat(n_, nlist_, arma::fill::zeros);
    for (std::size_t j = 0; j < nlist_; ++j) {
      const arma::mat A = A_list_.subcube(0,0,j,count_in_design_-1,count_in_design_-1,j);
      const arma::vec u_idx_in = u_list_.col(j);
#pragma omp parallel for
      for (std::size_t i = 0; i < n_; ++i) {
        val_out_mat(i,j) = remove_one_many(A, match_uvec(exp_cond_in_design_,idx_in_(i)), u_idx_in);
      }
    }
    
    return rd_mode_==1 ? val_out_mat * weights_ : arma::vec(arma::max(val_out_mat, 1));
  }
  
  
  double add_one_helper_many(const arma::mat &A, 
                             arma::uword i,
                             const arma::uvec &idx_orig, 
                             arma::uvec ii) {
    // ii is the rows of sigma for the new condition to add
    // i is the index of the design
    arma::uword n_to_add = ii.n_elem; // number of observations to add
    arma::uword n_already_in = idx_orig.n_elem; // number of observations currently in the design
    arma::uvec idx_in_vec(n_already_in + n_to_add, fill::zeros);
    idx_in_vec(span(0,n_already_in - 1)) = idx_orig;
    arma::mat A2 = A;
    
    // loop through all the indices to add
    for(arma::uword j=0; j<n_to_add; j++){
      A2 = add_one_mat(A2, 
                       sig_list_.slice(i)(ii(j),ii(j)),
                       sig_list_.slice(i).submat(idx_in_vec(span(0,n_already_in - 1)), arma::uvec({ii(j)})));
      idx_in_vec(n_already_in ) = ii(j);
      n_already_in += 1;
    }
    
    arma::uvec keep_idx = idx_in_vec(span(0,n_already_in - 1));
    arma::vec ucol = u_list_.col(i);
    return obj_fun(A2,ucol.elem(keep_idx));
  }
  
  arma::mat add_one_helper_many_mat(const arma::mat &A, 
                                    arma::uword i,
                                    const arma::uvec &idx_orig, 
                                    arma::uvec ii) {
    
    // ii is the rows of sigma for the new condition to add
    arma::uword n_to_add = ii.n_elem; // number of observations to add
    arma::uword n_already_in = idx_orig.n_elem; // number of observations currently in the design
    arma::uvec idx_in_vec(n_already_in + n_to_add, fill::zeros);
    idx_in_vec(span(0,n_already_in - 1)) = idx_orig;
    arma::mat A2 = A;
    
    // loop through all the indices to add
    for(arma::uword j=0; j<n_to_add; j++){
      
      A2 = add_one_mat(A2, 
                       sig_list_.slice(i)(ii(j),ii(j)),
                       sig_list_.slice(i).submat(idx_in_vec(span(0,n_already_in - 1)), arma::uvec({ii(j)})));
      
      idx_in_vec(n_already_in) = ii(j);
      
      n_already_in += 1;
    }
    return A2;
  }
  
  // evaluate add_one
  // same changes as for the val_out comp function
  arma::vec comp_val_in_vec(bool use_rm1A = false) {
    arma::mat val_in_mat(idx_out_.n_elem,nlist_,arma::fill::zeros);
    for (arma::uword j = 0; j < nlist_; ++j) {
      const arma::mat A = use_rm1A ? rm1A_list_.subcube(0,0,j,count_in_rm1_-1,count_in_rm1_-1,j) : 
      A_list_.subcube(0,0,j,count_in_design_-1,count_in_design_-1,j);
#pragma omp parallel for
      for (arma::uword i = 0; i < n_cond_ - n_; ++i) {
        const arma::uvec rows_keep = use_rm1A ? rows_in_rm1_(span(0,count_in_rm1_-1)) : 
        rows_in_design_(span(0,count_in_design_-1));
        val_in_mat(i,j) = add_one_helper_many(A, j, rows_keep, match_uvec(exp_cond_,idx_out_(i)));
      }
    }
    return rd_mode_==1 ? val_in_mat * weights_ : arma::vec(arma::max(val_in_mat, 1));
  }
  
  
  //REORDER FUNCTION
  void reorder_obs() {
    
    // find one index from idx_in to remove
    // which results in largest val of remove_one()
    arma::vec val_out_vec = comp_val_out_vec(); 
    idx_in_ = idx_in_.rows(arma::sort_index(val_out_vec, "descend"));
    //if(trace_)Rcpp::Rcout << "done comp out vec\n";
    // now update the indexes to ignore if they are fixed
    if(any_fix_ > 0){
      arma::uvec ignore(n_,fill::zeros);
      arma::uword count_ignore = 0;
      for(arma::uword j=0; j<n_; j++){
        bool any_to_ignore = any(nfix_ == idx_in_(j));
        if(any_to_ignore){
          ignore(count_ignore) = j;
          ++count_ignore;
        }
      }
      ignore_idx_ = ignore.rows(0,count_ignore-1);
    }
    // find one index from idx_out to add (swap)
    // which results in largest val of add_one()
    arma::vec val_in_vec = comp_val_in_vec();
    idx_out_ = idx_out_.rows(arma::sort_index(val_in_vec, "descend"));
    //if(trace_)Rcpp::Rcout << "done comp in vec\n";
  }
};



// [[Rcpp::export]]
Rcpp::List GradRobustStep(arma::uword N,
                          arma::uvec idx_in, 
                          Rcpp::List C_list, 
                          Rcpp::List X_list, 
                          Rcpp::List sig_list, 
                          arma::vec weights,
                          arma::uvec exp_cond,
                          arma::uvec nfix, 
                          arma::uword any_fix = 0,
                          arma::uword rd_mode = 1,
                          bool trace = true) {
  
  // need to map Rcpp list to the field here:
  arma::uword ndesign = weights.n_elem;
  arma::field<arma::vec> Cfield(ndesign);
  arma::field<arma::mat> Xfield(ndesign);
  arma::cube sigcube(N,N,ndesign);
  for(arma::uword j=0; j<ndesign; j++){
    Cfield[j] = as<arma::vec>(C_list[j]);
    Xfield[j] = as<arma::mat>(X_list[j]);
    sigcube.slice(j) = as<arma::mat>(sig_list[j]);
  }
  
  HillClimbing hc(idx_in, Cfield, Xfield, sigcube, weights, exp_cond, any_fix, nfix, rd_mode, trace);
  hc.grad_robust_step();
  return Rcpp::List::create(Named("idx_in") = hc.idx_in_,
                            Named("idx_out") = hc.idx_out_,
                            Named("best_val_vec") = hc.best_val_vec_);
}

// [[Rcpp::export]]
Rcpp::List GradRobustAlg1(arma::uword N,
                          arma::uvec idx_in, 
                          Rcpp::List C_list, 
                          Rcpp::List X_list, 
                          Rcpp::List sig_list, 
                          arma::vec weights,
                          arma::uvec exp_cond,
                          arma::uvec nfix, 
                          arma::uword any_fix = 0,
                          arma::uword rd_mode = 1,
                          bool trace = true) {
  
  // need to map Rcpp list to the field here:
  arma::uword ndesign = weights.n_elem;
  arma::field<arma::vec> Cfield(ndesign);
  arma::field<arma::mat> Xfield(ndesign);
  arma::cube sigcube(N,N,ndesign);
  for(arma::uword j=0; j<ndesign; j++){
    Cfield[j] = as<arma::vec>(C_list[j]);
    Xfield[j] = as<arma::mat>(X_list[j]);
    sigcube.slice(j) = as<arma::mat>(sig_list[j]);
  }
  
  HillClimbing hc(idx_in, Cfield, Xfield, sigcube, weights, exp_cond, any_fix, nfix, rd_mode, trace);
  hc.grad_robust_alg1();
  return Rcpp::List::create(Named("idx_in") = hc.idx_in_,
                            Named("idx_out") = hc.idx_out_,
                            Named("best_val_vec") = hc.best_val_vec_);
}
