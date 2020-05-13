
bool search(std::vector<int>* match_loc,std::vector<int>* num_match,std::string find,std::string file_address,int start_pos=0,int stop_pos=0);
int molp_sym_parser(std::string file);
int molp_method_parser(std::vector<int>* method_pos,std::string file);
bool molp_wf_parser(const int method_index,std::vector<int>* n_elec,std::vector<int>* sym,std::vector<int>* spin,std::vector<int>*charge,std::vector<int>* n_states,std::string file);
bool molp_cas_reader(int method_index,std::vector<int>* n_occ,std::vector<int>* n_closed,std::vector<int>* n_frozen,std::string file);
bool molp_basis_size_parser(std::vector<int>* basis_size,std::string file);
bool molp_basis_parser(std::vector<int>* basis_size,std::vector<int>* cont_num,std::vector<int>* nuc_bas_func,std::vector<int>* l_val,std::vector<int>* m_val,std::vector<double>* cont_zeta,std::vector<double>* cont_coeff,std::string file);
