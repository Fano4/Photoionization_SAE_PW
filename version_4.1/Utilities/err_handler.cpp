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
   std::cout<<"ERROR. MORE THAN ONE BASIS DATA FOUND IN INPUT FILE "<<std::endl<<file.c_str()<<std::endl;
   exit(EXIT_SUCCESS);
}
