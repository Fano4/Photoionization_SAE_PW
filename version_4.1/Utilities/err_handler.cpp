#include <cstdlib>
#include <iostream>
#include <string>

void err_file_not_found(std::string file)
{
   std::cout<<"ERROR. FILE NOT FOUND"<<std::endl<<file.c_str()<<" DOES NOT EXIST. EXIT"<<std::endl;
   exit(EXIT_SUCCESS);
}
void err_sym_not_found(std::string file)
{
   std::cout<<"ERROR. SYM CARD NOT FOUND IN INPUT FILE "<<std::endl<<file.c_str()<<std::endl;
   exit(EXIT_SUCCESS);
}
void err_sym_not_rec(std::string file,std::string symdata)
{
   std::cout<<"ERROR. SYMMETRY POINT GROUP NOT RECOGNIZED : "<<symdata.c_str()<<std::endl;
   std::cout<<" ENCOUTERED IN FILE "<<file.c_str()<<std::endl;

   exit(EXIT_SUCCESS);
}
void err_input_end_not_found(std::string file)
{
   std::cout<<"ERROR. END OF THE INPUT FILE NOT FOUND"<<std::endl<<file.c_str()<<std::endl;
   exit(EXIT_SUCCESS);
}
void err_wf_not_found(std::string file,int method_index)
{
   std::cout<<"ERROR. WF CARD NOT FOUND IN INPUT FILE "<<std::endl<<file.c_str()<<std::endl;
   std::cout<<"PARSING METHOD INDEX  "<<method_index<<std::endl;
   exit(EXIT_SUCCESS);
}
void err_wf_too_many_param(std::string file,std::string parsed)
{
   std::cout<<"ERROR. TOO MANY WF PARAMETERS IN INPUT FILE "<<std::endl<<file.c_str()<<std::endl;
   std::cout<<"WHILE PARSING THE WF FLAG  "<<parsed<<std::endl;
   exit(EXIT_SUCCESS);
}
void err_cas_too_many_ir(std::string file,std::string parsed)
{
   std::cout<<"ERROR. TOO MANY IRREDUCIBLE REPRESENTATIONS IN INPUT FILE "<<std::endl<<file.c_str()<<std::endl;
   std::cout<<"WHILE PARSING THE OCC/CLOSED/FROZEN FLAG  "<<parsed<<std::endl;
   exit(EXIT_SUCCESS);
}
void err_basis_not_found(std::string file)
{
   std::cout<<"ERROR. BASIS DATA NOT FOUND IN INPUT FILE "<<std::endl<<file.c_str()<<std::endl;
   exit(EXIT_SUCCESS);
}
void err_too_many_basis_data(std::string file)
{
   std::cout<<"WARNING MORE THAN ONE BASIS DATA FOUND IN INPUT FILE "<<std::endl<<file.c_str()<<std::endl;
   std::cout<<"READING THE FIRST INSTANCE OF BASIS SET"<<std::endl;
}
void err_lcao_too_many_method(std::string file)
{
   std::cout<<"ERROR MORE THAN ONE METHOD FOUND IN INPUT FILE "<<std::endl<<file.c_str()<<std::endl;
   std::cout<<"CANNOT PARSE MORE THAN ONE LCAO ARRAY"<<std::endl;
}
void err_lcao_method_not_supported(int method,int method_pos,std::string file)
{
   std::cout<<"ERROR METHOD NOT SUPPORTED FOR COMPUTING LCAO COEFF IN INPUT FILE "<<std::endl<<file.c_str()<<std::endl;
   std::cout<<"HIT THE ERROR WHEN PARSING METHOD "<<method<<" AT POSITION "<<method_pos<<std::endl;
   exit(EXIT_SUCCESS);
}
void err_lcao_not_found(std::string file)
{
   std::cout<<"ERROR. LCAO COEFFICIENT BLOCK NOT FOUND IN INPUT FILE "<<std::endl<<file.c_str()<<std::endl;
   exit(EXIT_SUCCESS);
}
void err_civector_not_found(std::string file)
{
   std::cout<<"ERROR. CI VECTOR BLOCK NOT FOUND IN INPUT FILE "<<std::endl<<file.c_str()<<std::endl;
   exit(EXIT_SUCCESS);
}
void err_ci_not_rec_symbol(std::string csf_str,const char symbol)
{
   std::cout<<"ERROR. CSF SYMBOL NOT RECOGNIZED : "<<symbol<<std::endl;
   std::cout<<" ENCOUTERED IN STRING "<<csf_str.c_str()<<std::endl;

   exit(EXIT_SUCCESS);
}
void err_geom_not_found(std::string file)
{
   std::cout<<"ERROR. ATOMIC COORDINATES BLOCK NOT FOUND IN INPUT FILE "<<std::endl<<file.c_str()<<std::endl;
   exit(EXIT_SUCCESS);
}
void err_too_many_geom(std::string file)
{
   std::cout<<"WARNING: MULTIPLE ATOMIC COORDINATES BLOCKS IN INPUT FILE "<<std::endl<<file.c_str()<<std::endl;
   std::cout<<"READING 1ST ATOMIC COORDINATES BLOCK"<<std::endl;
}
