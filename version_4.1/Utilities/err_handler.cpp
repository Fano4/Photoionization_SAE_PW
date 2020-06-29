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
void err_lcao_not_found(std::string file,std::string search)
{
   std::cout<<"ERROR. LCAO COEFFICIENT BLOCK NOT FOUND IN INPUT FILE "<<std::endl<<file.c_str()<<std::endl;
   std::cout<<" SEEKING FOR STRING "<<search<<std::endl;
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
void err_diff_on_vec_size(int size1,int size2)
{
   std::cout<<"ERROR. OVERLAP CAN NOT BE COMPUTED BETWEEN CONFIGURATIONS OF DIFFERENT DIMENSIONS"<<std::endl;
   std::cout<<size1<< " != "<<size2<<std::endl;
   exit(EXIT_SUCCESS);
}
void err_end_of_file(std::string file,std::string method)
{
   std::cout<<"ERROR. UNEXPECTEDLY REACHED THE END OF INPUT FILE "<<file.c_str()<<std::endl;
   std::cout<<"ERROR RAISED WHILE IN "<<method.c_str()<<std::endl;
   exit(EXIT_SUCCESS);
}
void err_caution_after_geom(std::string file,int pos)
{
   std::cout<<"WARNING: CAUTION MESSAGE IN MOLPRO OUTPUT FILE AFTER GEOMETRY BLOCK. RISK OF UNDEFINED BEHAVIOUR IF THE MESSAGE IS SERIOUS"<<std::endl;
   std::cout<<"ERROR RAISED WHILE IN "<<file.c_str()<<" AT POSITION "<<pos<<std::endl;
}
void err_bad_indices_gen_I_integ(int l1,int l2)
{
   std::cout<<"ERROR: BAD INDICES USED FOR COMPUTING GENERALIZED I INTEGRAL. EXITING"<<std::endl;
   std::cout<<"THE INDICES ARE L1 = "<<l1<<" AND L2 = "<<l2<<std::endl;
   exit(EXIT_SUCCESS);
}
