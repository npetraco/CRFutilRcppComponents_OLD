phi_vec.ones();

// ff test code
arma::Col<int> f10 = ff_C(1);
arma::Col<int> f01 = ff_C(2);

Rcout << f10 << endl;
Rcout << f01 << endl;




// old node par fix code
arma::Row<int> anode_par_vec(2);
IntegerVector tmpn(2);
arma::Mat<int> node_par_new();

// edge par fix old
arma::Mat<int> aedge_par_mat(2,2);
IntegerVector tmpe(4);
List edge_par_new(edge_par.size());
for(int i=0; i<edge_par.size(); ++i) {
  
  tmpe = as<IntegerVector>(edge_par(i));
  aedge_par_mat(0,0) = tmpe(0);
  aedge_par_mat(0,1) = tmpe(1);
  aedge_par_mat(1,0) = tmpe(2);
  aedge_par_mat(1,1) = tmpe(3);
  
  edge_par_new(i) = aedge_par_mat;
  
  //Rcout << aedge_par_mat << endl;
  
}



// // Remove third index from node_par
// arma::Mat<int> node_par_new(node_par.n_rows,2);
// node_par_new = node_par.slice(0);
// Rcout << node_par_new << endl;
// 
// // Remove third index from each of the edge_par
// arma::Mat<int> aedge_par_mat(2,2);
// IntegerVector tmpe(4);
// List edge_par_new(edge_par.size());
// for(int i=0; i<edge_par.size(); ++i) {
//   
//   tmpe = as<IntegerVector>(edge_par(i));
//   aedge_par_mat(0,0) = tmpe(0);
//   aedge_par_mat(0,1) = tmpe(1);
//   aedge_par_mat(1,0) = tmpe(2);
//   aedge_par_mat(1,1) = tmpe(3);
//   
//   edge_par_new(i) = aedge_par_mat;
//   
//   Rcout << aedge_par_mat << endl;
//   
// }
// 
// List theta_pars;
// theta_pars["node_par"] = node_par_new;
// theta_pars["edge_par"] = edge_par_new;


// REMOVE THIS EVENTUALLY TO A SEPARATE FUNCTION SO WE ONLY HAVE TO DO IT ONCE IF THIS IS CALLED MANY TIMES
// EVENTUALLY WE ONLY WANT TO SEND IN A REF TO THESE IN THE "2D" FORMAT......... 
//List theta_pars         = fix_node_and_edge_par(node_par_crf,edge_par_crf);
//arma::Mat<int> node_par = theta_pars["node_par"];
//List           edge_par = theta_pars["edge_par"];



Rcout << "Edge: " << i << " "<< edge_mat(i,0) << "---" << edge_mat(i,1) << endl;
Rcout << "Left spin: "  << config(left_off) << ", Right spin: " << config(right_off) << endl;
Rcout << endl;
Rcout << ff_C(config(left_off)) << aepm << ff_C(config(right_off)) << endl;
Rcout << endl;
Rcout << "tmp: " << tmp << "phi_off: " << phi_off << endl;
Rcout << "===========================" << endl;
