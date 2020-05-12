#include <cstdlib>
#include <iostream>
#include <string>
#include <sstream>
#include <fstream>
#include <vector>
#include "utilities.h"

bool search(std::vector<int>* match_loc,std::vector<int>* num_match,std::string find,std::string file_address,int start_pos=0,int stop_pos=0)
{
   //The search algorithm open the file and put the buffer into a tring. The keyword is searched inside the string until eof.
   //Upon match, the location of the first charater of the matching result is stored in a vector
   //Upon eof, the vector is translated into an array and its length is recorded.
   //Return if 0 if not found, and 1 if found

    using namespace std;
    int count(0);
    size_t found;
    std::vector<int> match_loc_vec;

    ostringstream tmp_str;
    string str;
    ifstream file;
    
    //extract the text and put it in a string
    file.open(file_address.c_str());
    if (!file.is_open())
    {
        std::cout<<"Issue while opening text file at address "<<file_address.c_str()<<endl;
        exit(EXIT_SUCCESS);
    }
    else
    {
       tmp_str<<file.rdbuf();
       str=tmp_str.str();
    }
    file.close();

    // If a stop position is set, then extract a truncated substrin
    if(stop_pos!=0)
       str=str.substr(0,stop_pos);
    //search for the string find in the string str. We start by one seek, and then we loop while the result is not found.
    found = str.find(find,start_pos); 

    //if not found, return 0
    if (found == string::npos) 
        return 0;

    //if it was found, carry on and count until there is none left
    else
    {
      count++;
      match_loc->push_back(found);
       while(found!=string::npos)
       {
          found = str.find(find,found+1);
          if(found!=string::npos)
          {
             match_loc->push_back(found);
             count++;
          }

       }
       num_match->push_back(count);
       return 1;
    }
}


int molp_sym_parser(std::string file)
{

   //This function parses the symmetry statement of the molpro input file and returns the number of irreducible representations
   using namespace std;

   vector<int> position;
   vector<int> stop_pos;
   vector<int> match_num;

   string line;
   string data_str;

   ifstream input;

   //Search the end of the input section
   if(!search(&stop_pos,&match_num,"Passed",file))
      err_input_end_not_found(file);
   //Search in the file for the different symmetry keywords
   if(search(&position,&match_num,"symmetry",file,0,stop_pos[0]))
   {
      //If we match for symmetry, then we extract the first matching line 
      input.open(file.c_str());
      if(input.is_open())
      {
         input.seekg(position[0]);
         getline(input,line);
         input.close();
      }
      else 
         err_file_not_found(file);

      //Once we've got the line, we can parse it. We extract the data_substring by removing the symmetry keyword and the first comma

      data_str=line.substr(9); //remove "symmetry," 
   }
   else if(search(&position,&match_num,"sym",file,0,stop_pos[0]))
   {
      //If we match for sym, then we extract the first matching line 
      input.open(file.c_str());
      if(input.is_open())
      {
         input.seekg(position[0]);
         getline(input,line);
         input.close();
      }
      else 
         err_file_not_found(file);

      //Once we've got the line, we can parse it. We extract the data_substring by removing the symmetry keyword and the first comma

      data_str=line.substr(4); //remove "sym," 
   }
   else
      err_sym_not_found(file);

   //Seek for nosym keyword. This correspond to 1 IR
   if(data_str.find("nosym") != string::npos)
      return 1;


   //Seek for x,y or z keyword. This returns 2 IR for Cs
   else if( (data_str.find("x") != string::npos) ^ (data_str.find("y") != string::npos) ^ (data_str.find("z") != string::npos))
      return 2;
      
   //Seek for xy,xz or yz keyword. This returns 2 IR for C2 
   else if((data_str.find("xy") != string::npos) ^ (data_str.find("xz") != string::npos) ^ (data_str.find("yz") != string::npos))
      return 2;

   //Seek for xyz keyword. This returns 2 IR for Ci 
   else if(data_str.find("xyz") != string::npos)
      return 2;


   //Seek for "x,y","y,z","x,z" keyword. This returns 4 IR for C2v
   else if((data_str.find("x,y") != string::npos) ^ (data_str.find("x,z") != string::npos) ^ (data_str.find("y,z") != string::npos))
      return 4;

   //Seek for "xy,z","yz,x","xz,y" keyword. This returns 4 IR for C2h
   else if((data_str.find("xy,z") != string::npos) ^ (data_str.find("yz,x") != string::npos) ^ (data_str.find("xz,y") != string::npos))
      return 4;

   //Seek for "xz,yz","xy,yz","xy,xz" keyword. This returns 4 IR for D2
   else if((data_str.find("xy,yz") != string::npos) ^ (data_str.find("xz,yz") != string::npos) ^ (data_str.find("xy,xz") != string::npos))
      return 4;

   //Seek for x,y,z keyword. This returns 8 IR for D2h
   else if(data_str.find("x,y,z") != string::npos)
      return 8;

   else
      err_sym_not_rec(file,data_str);
   
   return 0; // Dead end
}




int molp_method_parser(std::vector<int>* method_pos,std::string file)
{
      //search for the different methods that are used in the output file
      //for now, only casscf/mcscf is supported. it will be possible to support other methods in the future.
      //The methods that are found are registered with the position of their block in the file. 
      //This will be used in the following functions

   std::vector<int> position;
   std::vector<int> stop_pos;
   std::vector<int> match_num;
   std::vector<int> method;
   int len;

   //Search the end of the input section
   if(!search(&stop_pos,&match_num,"Passed",file))
      err_input_end_not_found( file);
   // rhf = 1
   if(search(&position,&match_num,"rhf",file,0,stop_pos[0])) 
   {
      for(int i=0;i!=match_num[0];i++)
         method.push_back(1);
   }
   match_num.clear();
   //casscf = 2
   if(search(&position,&match_num,"casscf",file,0,stop_pos[0])) 
   {
      for(int i=0;i!=match_num[0];i++)
         method.push_back(2);
   }
   match_num.clear();
   //mrci = 3
   if(search(&position,&match_num,"mrci",file,0,stop_pos[0])) 
   {
      for(int i=0;i!=match_num[0];i++)
         method.push_back(3);
   }
   match_num.clear();
   len=position.size();

   for(int i=0;i!=len;i++)
   {
       method_pos->push_back(method[i]);
       method_pos->push_back(position[i]);
   }

   return 0;
}





bool molp_wf_parser(const int method_index,std::vector<int>* n_elec,std::vector<int>* sym,std::vector<int>* spin,std::vector<int>*charge,std::vector<int>* n_states,std::string file)
{
   //Here, we parse the wf commands of the program blocks. The function allows to chose the index of the method that has been parsed.
   //Symmetry must be given as an input since it conditions the number of wf statements we can meet in one command block
   //The vectors will be pushed back with the elements associated with all the irreducible representations for the command block
   //We do not know the size of the vectors in advance since the number of wf cards depend on the number of spin states
   //If the wf statement for a given IR is not set, its n_states is set to zero

   using namespace std;

   vector<int> method_pos;
   vector<string> wf_cards;
   ifstream input;
   string line;
   size_t pos;


   //Recall metho parser to get the actual position of the command block we want to parse
   molp_method_parser(&method_pos,file);

   //Once we are in the command block, we parse the block line by line until we hit a closing brace.
   input.open(file.c_str());
   if(!input.is_open())
      err_file_not_found(file);

   // Get to the first line of the command block
   input.seekg(method_pos[2*method_index+1]); 

   bool eob(0); //end_of_block
   while(!eob)
   {
      //extract line by line
      getline(input,line);

      //save line whenever a wf instance is found,excluding the "wf," card.
      pos=line.find("wf");
      if(pos!=string::npos)
         wf_cards.push_back(line.substr(3));

      //check for the closing brace symbol
      pos=line.find("}");
      if(pos!=string::npos)
         eob=1;
   }
   input.close();

   //Check is at least one wf card has been found
   if( wf_cards.size() == 0 )
      err_wf_not_found(file,method_index);

   input.close();

   string temp;
   size_t pos2;
   std::vector<int> pos_commas;
   stringstream ss;
   unsigned int count;
   int temp_int;
   bool test(0);
   bool test1(0);

   //Now we read the content of the wf cards. There should be one per symmetry and spin.
   for(unsigned int i=0;i!= wf_cards.size();i++)
   {
      //std::cout<<i<<std::endl;
      ss.str(wf_cards[i]);
      getline(ss,temp);

      //seek for state card
      pos=wf_cards[i].find(";state,");

      //if state card found, save the number of states and remove the tail of the wf card
      if(pos!=string::npos)
      {
         ss.str(wf_cards[i].substr(pos+7));
         pos=ss.str().find(';');
         if(pos==string::npos)
         {
            std::cout<<"Warning: no ';' found after number of states specifications. Risk of unexpected behaviour"<<std::endl;
            n_states->push_back(0);
            break;
         }
         else
            temp=ss.str().substr(pos-2,pos);
         temp_int=atoi(temp.c_str());
         n_states->push_back(temp_int);
         pos=wf_cards[i].find(";state,");
         temp=wf_cards[i].substr(0,pos);
         wf_cards[i]=temp;
      }

      //Browse the content of the wf vector
      ss.str(wf_cards[i]);
      
      // Numbers can range from 0 to 999
       test=0;
       test1=0;
      count=0;
      while(!test)
      {
      //   std::cout<<temp<<endl;
         pos=ss.str().find(',');
         pos2=ss.str().find(',',pos+1);
         if(pos2==string::npos)
         {
            temp=ss.str().substr(pos+1);
            test1=1;
         }
         else
            temp=ss.str().substr(pos+1,pos2-pos);


         //Each number correspond to a definite quantity. 
         // 0 : n_elec , 1 : sym, 2 : spin , 3 : charge
         switch(count)
         {
            case 0:
               n_elec->push_back(atoi(temp.c_str()));
               break;
            case 1:
               sym->push_back(atoi(temp.c_str()));
               break;
            case 2:
               spin->push_back(atoi(temp.c_str()));
               if(test1==1)
                  test=1;
               break;
            case 3:
               charge->push_back(atoi(temp.c_str()));
               if(test1==1)
                  test=1;
               break;

            default:
               err_wf_too_many_param(file,wf_cards[i]);
               break;
         }
         if(pos2==string::npos)
            temp="";
         else
            temp=ss.str().substr(pos2);

         ss.str(temp);

         count++;
      }

   }

   return 0;
}
bool molp_cas_reader(int method_index,std::vector<int>* n_occ,std::vector<int>* n_closed,std::vector<int>* n_frozen,std::string file)
{
   //Extracting the active space parameters in a given command block

   using namespace std;

   vector<int> method_pos;
   ifstream input;
   string line;
   string occ_str;
   string closed_str;
   string frozen_str;
   size_t pos;
   bool occ_check,closed_check,frozen_check;
   unsigned int n_sym;


   //Recall metho parser to get the actual position of the command block we want to parse
   molp_method_parser(&method_pos,file);
   n_sym=molp_sym_parser(file);
   bool eob(0); //end_of_block

   input.open(file.c_str());
   if(!input.is_open())
      err_file_not_found(file);

   // Get to the first line of the command block
   input.seekg(method_pos[2*method_index+1]); 
   while(!eob)
   {

      //extract line by line
      getline(input,line);

      //save value to n_occ upon n_occ match
      pos=line.find("occ");
      if(pos!=string::npos)
      {
         occ_str=line.substr(4).c_str();
         occ_check=1;
      }
      
      //save value to n_closed upon n_closed match
      pos=line.find("closed");
      if(pos!=string::npos)
      {
         closed_str=line.substr(7).c_str();
         closed_check=1;
      }

      //save value to n_frozen upon n_frozen match
      pos=line.find("frozen");
      if(pos!=string::npos)
      {
         frozen_str=line.substr(7).c_str();
         frozen_check=1;
      }

//      std::cout<<line<<std::endl;
      //check for the closing brace symbol
      pos=line.find("}");
      if(pos!=string::npos)
         eob=1;
   }
      //std::cout<<occ_str<<","<<closed_str<<","<<frozen_str<<std::endl;

   if(!occ_check)
         std::cout<<" Warning: n_occ not set in active space "<<std::endl<<"method index "<<method_index<<std::endl<<"file "<<file.c_str()<<endl;

   if(!closed_check)
         std::cout<<" Warning: n_closed not set in active space "<<std::endl<<"method index "<<method_index<<std::endl<<"file "<<file.c_str()<<endl;
   input.close();

   //Parse the strings for occ, closed and frozen. They should have a length = n_sym
   stringstream ss;
   string temp;
   size_t pos2;
   unsigned int count(0);

   if(occ_check)
   {
      count=0;
         ss.str(occ_str);
      while( count < n_sym )
      {
         pos=ss.str().find(',');
         pos2=ss.str().find(',',pos+1);
         if(pos2==string::npos)
            temp=ss.str().substr(pos+1,pos+2);
         else
         {
            temp=ss.str().substr(pos+1,pos2-pos);
            ss.str(ss.str().substr(pos2));
         }

         if(count<n_sym)
            n_occ->push_back(atoi(temp.c_str()));
         else
            err_cas_too_many_ir(file,temp);
         count++;
      }

   }
   if(closed_check)
   {
      count=0;
      ss.str(closed_str);
      while( count < n_sym )
      {
         pos=ss.str().find(',');
         pos2=ss.str().find(',',pos+1);
         if(pos2==string::npos)
            temp=ss.str().substr(pos+1,pos+2);
         else
         {
            temp=ss.str().substr(pos+1,pos2-pos);
            ss.str(ss.str().substr(pos2));
         }

         if(count<n_sym)
            n_closed->push_back(atoi(temp.c_str()));
         else
            err_cas_too_many_ir(file,temp);
         count++;
      }

   }
   if(frozen_check)
   {
      count=0;
      ss.str(frozen_str);
      while( count < n_sym )
      {
         pos=ss.str().find(',');
         pos2=ss.str().find(',',pos+1);
         if(pos2==string::npos)
            temp=ss.str().substr(pos+1,pos+2);
         else
         {
            temp=ss.str().substr(pos+1,pos2-pos);
            ss.str(ss.str().substr(pos2));
         }

         if(count<n_sym)
            n_frozen->push_back(atoi(temp.c_str()));
         else
            err_cas_too_many_ir(file,temp);
         count++;
      }

   }
   

   return 0;
}
bool molp_basis_size_parser(int* basis_size,std::string file)
{
   using namespace std;

   int n_sym;
   size_t pos;
   ifstream input;
   string tmp_str;
   string line;

   n_sym=molp_sym_parser(file);

   input.open(file.c_str());
   if(!input.is_open())
      err_file_not_found(file);

   while(!input.eof())
   {
      //extract line by line
      getline(input,line);

      //save value to basis_size upon match
      pos=line.find("NUMBER OF CONTRACTIONS:");
      if(pos!=string::npos)
      {

         input.seekg(pos+23);
         input>>basis_size[0];

         if(n_sym>1)
         {
            input>>tmp_str;

            //If there is additional symmetry, save the size of the different IR to the basis size array
            for(int i=0;i!=n_sym;i++)
            {
               input>>tmp_str;
               basis_size[i]=atoi(separateThem(tmp_str).c_str());
               std::cout<<"sym "<<i+1<<" = "<<basis_size[i]<<std::endl;
               input>>tmp_str;
            }

         }
      }
      
   }
   input.close();
   return 0;
}
bool molp_basis_parser(int* basis_size,std::vector<int>* cont_num,std::vector<int>* nuc_bas_func,std::vector<int>* l_val,std::vector<int>* m_val,std::vector<double>* cont_zeta,std::vector<double>* cont_coeff,std::string file)
{
   using namespace std;


   ifstream input;



   // Find the basis data part and read the number of contractions for each basis function. 
   // The molpro file sorts the basis functions by symmetry and gives the contraction coefficients and zeta for all the basis fucntions.
   // Only support Pople type basis sets.

   vector<int> basis_position;
   vector<int> num_of_match;
   unsigned int count;
   int cur_line_pos;
   int prev_line_pos;
   int bas_func_index;
   int n_sym;
   size_t pos;
   string teststring;
   string tmp_str;
   string line;

   n_sym=molp_sym_parser(file);
   molp_basis_size_parser(basis_size,file);
   // Search for the basis data block in the input file
   search(&basis_position,&num_of_match,"BASIS DATA",file);
   stringstream ss,sstream;

   // Check if the basis data block was found and if it is unique in the file
   if( num_of_match.size() < 1 )
      err_basis_not_found(file);
   else if(num_of_match[0] > 1 )
      err_too_many_basis_data(file);

   //If the basis data is OK, then proceed reading the basis set data.
   else
   {
      input.open(file);
      if(!input.is_open())
         err_file_not_found(file);

      //Start where the BASIS DATAS entry matches and seek for the line where the first basis function is defined.
      input.seekg(basis_position[0]);
      cur_line_pos=input.tellg();
      do
      {
         getline(input,line);
         pos=line.find("1.1");
         //keep the position of the lines so that we can go back after reading the matching file
         prev_line_pos=cur_line_pos;
         cur_line_pos=input.tellg();
      }while(pos == string::npos);

      //Go back to the first line of basis set definition
      input.seekg(prev_line_pos);

      //The basis function belong to the different IR of the point group. We browse all of them by symmetry
      bas_func_index=0;
      for(int s=0;s!=n_sym;s++)
      {
         //for each IR, we know the number of basis functions 
         for(int i=0;i!=basis_size[s];i++)
         {

            //We read one basis function at a time. We don't know how many contractions there are for a single basis func. 
            //We stop when we reach the next basis function, whic is defined below
            sstream.str("");
            if(i<basis_size[s]-1)
            {
               sstream<<i+2<<"."<<s+1;
               teststring=sstream.str();
            }

            //If we reach the end of an IR block, we move on with the first ao on the next block
            else if(i == basis_size[s]-1 && s<n_sym-1)
            {
               sstream<<"1."<<s+2;
               teststring=sstream.str();
            }

            //if we reach the last element of the last block, we stop at the next statement which is the keyword "NUCLEAR"
            else if (i==basis_size[s]-1 && s == n_sym-1)
               teststring="NUCLEAR";

            //basis function label
            input>>tmp_str;
            //irreducible rep in ascii 
            input>>tmp_str;
            //nuc_bas_func
            input>>tmp_str;
            nuc_bas_func->push_back(atoi(tmp_str.c_str()));
            //angular numbers
            input>>tmp_str;

            l_val->push_back(l_number(tmp_str));
//            m_val->push_back(ml_number(tmp_str,*l_val[bas_func_index]));
            std::cout<<l_val->size()<<std::endl;
            m_val->push_back(ml_number(tmp_str,0));


            count=0;
            while(tmp_str!=teststring)
            {
               input>>tmp_str;
               if( !bool(count%2) && tmp_str!=teststring )
               {
                  cont_zeta->push_back(atof(tmp_str.c_str()));
                  count++;
               }
               else if(tmp_str!=teststring)
                  cont_coeff->push_back(atof(tmp_str.c_str()));

               if(input.eof())
               {
                  break;
               }
            }
            cont_num->push_back(count);
            if(input.eof())
            {
               std::cout<<"UNEXPECTED EOF WHILE PARSING BASIS SET. EXIT"<<std::endl;
               exit(EXIT_SUCCESS);
            }
            bas_func_index++;
         }
      }
      input.close();
   }
   
   return 0;
}
