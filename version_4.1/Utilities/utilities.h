//err_handler
void err_file_not_found(std::string file);
void err_sym_not_rec(std::string file,std::string symdata);
void err_sym_not_found(std::string file);
void err_wf_not_found(std::string file,int method_index);
void err_wf_too_many_param(std::string file,std::string parsed);
void err_cas_too_many_ir(std::string file,std::string parsed);
void err_basis_not_found(std::string file);
void err_too_many_basis_data(std::string file);
void err_input_end_not_found(std::string file);
void err_lcao_method_not_supported(int method,int method_pos,std::string file);
void err_lcao_too_many_method(std::string file);
void err_lcao_not_found(std::string file);
//string_utils
enum charTypeT{ other, alpha, digit};
std::string separateThem(std::string inString);
charTypeT charType(char c);

int l_number(std::string bas_func_type);
int ml_number(std::string bas_func_type,int l);
