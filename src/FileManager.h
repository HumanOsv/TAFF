
#include "define.h"
#include <iomanip>
#include <sys/types.h>
#include <signal.h>
#include <limits.h>

using namespace std;

class fileManager {
	struct lines{				//list with info for basin, integral, volumne, XYZ coordinates and value.
		string basin, integral, volumen, c_x, c_y, c_z, value,att;
		lines *next;
	};
	struct deg_list{					//lst with degenerate data.
		string deg_line;
		deg_list *next;
	}	;
	lines *root, *root3;
	deg_list *root2;
	string sum, Multiwfn_Route;
	public:
		fileManager();//constructor de la clase...
		int Multiwfn_Pipeline(string, string, string);
		void RepresentationGraphics( string, int , string,int);
		void WriteMol2(string, int, string);
		void Extrac_N_Write_Integral_DAT(string, string, string);
		void TableWrite(string, int, int, int);
		string table(string , int );
		void DeleteFile(string);
		void MultiwfnRoute();
		string currentDateTime();
		void Density_Data_File(string, string, string);
		void PrintPcenOW(long double Occ[], int ,long double Unocc[], int, long double, long double, string );
};//fin definicion de la clase...
	
fileManager::fileManager(){
	root= new lines;
	root2= new deg_list;
	root3=new lines;	
}
void fileManager::Density_Data_File(string filename, string save, string aprox){
	ifstream file;
	file.open(filename.c_str());
	string data,palabra;
	int flag1=0,flag2=0, flag3=0;	//many flags to get specific data.
	lines *lector;
	root3->next=NULL;
	lector=root3;
 	int aux=0, aux2=0, radical_type=0;
	while(file){
		getline(file, data);
		if((strstr(data.c_str(),"X,Y,Z")!=NULL && strstr(data.c_str(),"Angstrom")!=NULL) || flag1==1){		//start found or already found it.
			if(strstr(data.c_str(),"grid")==NULL){				//end not found, so we wirte data in list . (basin, integral, volumen.)
				if(radical_type==1){			//if not degenerated attractor are found then we need to find the data a second and final time.
					aux=0;						//so we reset the aux and the flags
					flag2=0;
					radical_type=0;
				}
				if(flag1==1){
					istringstream mov(data);
					mov >> lector->basin >> lector->c_x>>lector->c_y>>lector->c_z >> lector->value;
					lector->next=new lines;
					lector=lector->next;
					aux++;
				}
				flag1=1;
			}else{								
				lector=root3;										
				flag1=0;
				flag2=1;
				radical_type=1;
			}	
		}
		if((strstr(data.c_str(),"Integral(a.u.)")!=NULL || flag3==1) && (flag2==1)){	//the first block of info has been written and the start of the second blocks has begun.
			if(strstr(data.c_str(),"Sum of above values")==NULL){						//the end isn't near.
				if(flag3==1 && aux2<aux){							
					istringstream iss(data);				//info writing.
					iss >>lector->basin >>lector->integral >>lector->volumen;
					lector=lector->next;
					aux2++;
				}
				flag3=1;
			}else{
				flag3=0;
				stringstream iss(data.substr(data.find(":")+1));	//1. we save the Sum of above values info.
				iss >> sum;											//2. we reset the list to write XYZ coordinates and value.				
			}
		}
	}
	file.close();		//close the tmp file.
	TableWrite(save, aux, 0	,2);
	RepresentationGraphics(save, aux, aprox,2);
}



void fileManager::DeleteFile(string file){
	remove(file.c_str());

}
int fileManager::Multiwfn_Pipeline(string orden, string file_to_open, string file_to_save){
	int state=0;
    ofstream output(file_to_save.c_str(),ios::app);
    // pipes for parent to write and read
    pipe(pipes[PARENT_READ_PIPE]);
    pipe(pipes[PARENT_WRITE_PIPE]);

    pid_t pid;
    pid=fork();
     
    if(pid==0) {
 
        dup2(CHILD_READ_FD, STDIN_FILENO);	//input
        dup2(CHILD_WRITE_FD, STDOUT_FILENO);	//output
        dup2(CHILD_WRITE_FD,STDERR_FILENO);		//error, to the same pipe of output.
 
        /* Close fds not required by child. Also, we don't
           want the exec'ed program to know these existed */
        close(CHILD_READ_FD);
        close(CHILD_WRITE_FD);
        close(PARENT_READ_FD);
        close(PARENT_WRITE_FD);
        //execl(Multiwfn_Route.c_str(),"./Multiwfn", file_to_open.c_str(), NULL);  //route, call MUltiwfn, order, NULL
        execlp("Multiwfn","Multiwfn", file_to_open.c_str(), NULL);  //route, call MUltiwfn, order, NULL
    } else {
        char buffer[100];
        int count;
 
        /* close fds not required by parent */       
        close(CHILD_READ_FD);
        close(CHILD_WRITE_FD);
 
        // Write to child’s stdin
        write(PARENT_WRITE_FD, orden.c_str(), orden.length());
        // Read from child’s stdout
        while((count=read(PARENT_READ_FD, buffer, sizeof(buffer)-1))){
            buffer[count] = 0;
            if(strstr(buffer,"attractor")!=NULL && strstr(buffer,"absolute value")!=NULL && strstr(buffer,"insignificant")!=NULL && orden=="17\n1\n2\nx\n"){
            	// Note: There has been a grid data in the memory, please select generating the basins by which manne
            	// Note: There are attractors having very low absolute value (<1.00E-05) and thus insignificant, how to deal with them?
            	state=1;
            	kill(pid, SIGKILL);
            }
            output << buffer ;	//write tmp files.
        }
    }
    output.close();
    return state;
}
void fileManager::Extrac_N_Write_Integral_DAT(string filename, string save, string aprox){
	ifstream file;
	file.open(filename.c_str());
	string data,palabra;
	int flag1=0,flag2=0, flag3=0, flag4=0, flag5=0;	//many flags to get specific data.
	lines  *lector;
	deg_list 	* deg_lec;
	root->next=NULL;
	lector=root;

	root2->next=NULL;
 	deg_lec=root2;
 	int aux=0, aux2=0, aux_d=0;
	while(file){
		getline(file, data);
		size_t f1=data.find("Integral(a.u.)");			//all the secrets passwords.
		size_t f2=data.find("Sum of above values:");
		size_t f3=data.find("Angstrom");
		size_t f4=data.find("The members");
		size_t f5=data.find("Basin analysis");
		if(f1 !=string::npos || flag1==1){		//start found or already found it.
			if(f2==string::npos){				//end not found, so we wirte data in list . (basin, integral, volumen.)
				if(flag1==1){//The members of degenerated attractor
					
					istringstream mov(data);
					mov >> lector->basin >> lector->integral>>lector->volumen;
					lector->next=new lines;
					lector=lector->next;
					aux++;
				}
				flag1=1;
			}else{								//if we found thw first blocks end then:	
				stringstream iss(data.substr(data.find(":")+1));	//1. we save the Sum of above values info.
				iss >> sum;											//2. we reset the list to write XYZ coordinates and value.
				lector=root;										//3. flags changes
				flag1=0;
				flag2=1;
			}	
		}
		if((f3 !=string::npos || flag3==1) && (flag2==1)){	//the first block of info has been written and the start of the second blocks has begun.
			if(f4==string::npos && f5==string::npos){						//the end isn't near.
				if(flag3==1 && aux2<aux){							
					istringstream iss(data);				//info writing.
					iss >> lector->att >>lector->c_x>>lector->c_y >>lector->c_z >> lector->value;
					lector=lector->next;
					aux2++;
				}
				flag3=1;
			}else{
				flag3=0;
				if(f4!=string::npos){
					flag4=1;
				}else{
					flag5=1;
				}
				
			}
		}
		if((f4 != string::npos || flag4==1) && flag5==0){			//last block of data (degenerate).
			if(f5==string::npos){
				deg_lec->deg_line=data;				//second list is fill.
				deg_lec->next=new deg_list;
				deg_lec=deg_lec->next;
				aux_d++;
			}else{
				flag4=0;
				flag5=1;
			}
		}
	}
	file.close();		//close the tmp file.
	TableWrite(save, aux, aux_d,1);
	RepresentationGraphics(save, aux, aprox,1);
	//WriteMol2(save,aux,aprox);
	
}

void fileManager::TableWrite(string save, int aux, int aux_d, int type){
	lines *lector;	
	deg_list *deg_lec;
	if(type==1){
		lector=root;
		save=save+"-integral.txt";
	}else{
		lector=root3;
		save=save+"-rho.txt";
	}

	ofstream output(save.c_str(),ios::trunc);	//output file, name may change.
	
	output << "\n" ;
	output << "\t  __________    ____    ________   ________\n" ;
	output << "\t  \\ __   ___\\  / __ \\  |  _____/  |  _____/\n" ;
	output << "\t      | |     / /  \\ \\ | |        | |\n"       ;	
	output << "\t      | |     | |__| | | |___     | |___\n"    ;	
	output << "\t      | |     |  __  | |  __/     |  __/\n"    ;	
	output << "\t      | |     | |  | | | |        | |\n"       ;	
	output << "\t      |_|     |_|  |_| |_|        |_|\n"       ;	
	output <<  "" ;	
	output << "\t# # # # # # # # # # # # # # # # # # # # # # # # # #\n" ;	
	output << "\t# Topological Analysis of Fukui Functions (TAFF)  #\n" ;	  
	output << "\t# DATE: "<<currentDateTime()<<"\t\t\t\t          #\n";	
	output << "\t# TiznadoLab\t\t\t\t\t\t\t\t\t  #\n";
	output << "\t# # # # # # # # # # # # # # # # # # # # # # # # # #\n" ;	 
	output <<  "\n" ;	
	output <<  "\n" ;	
	output << " ·······································································································\n" ;	
	output << " BASIN   VOLUME        INTEGRAL          <X>        <Y>         <Z>         VALUE\n" ;	
	output << " ·······································································································\n" ;	
	int i=0;
	cout.precision(7);
	if ( lector != 0 ) { //Makes sure there is a place to star
		while(i<aux){
		//while(lector->next!=NULL){
			i++;
			//output <<lector->basin << "\t"<< lector->volumen<<setw(3)<<setfill('0')<< "\t" << lector->integral<< "\t" << lector->c_x<< "\t" << lector->c_y<< "\t" << lector->c_z << "\t"<< lector->value << endl;
			
			output <<lector->basin << "\t"<< table(lector->volumen,12)<< "\t" << lector->integral<< "\t" << lector->c_x<< "\t" << lector->c_y<< "\t" << lector->c_z << "\t"<< lector->value << endl;
	    	
	    	lector = lector->next;
	  	}
	}
	output << " ·······································································································\n";	
	output << " TOTAL\t\t\t"<<sum<<"\n";	
	output << " \n" ;
	output << " \n" ;
	if(type==1){
		output << " ·······································································································\n" ;	
		output << " THE MEMBERS OF DEGENERATED ATTRACTOR\n";	
		output << " ·······································································································\n" ;	
		i=0;
		deg_lec=root2;
		if ( deg_lec != 0 ) { //Makes sure there is a place to start
			while(i<aux_d){
			//while(deg_lec->next!=NULL){
	  			i++;
	  			output << deg_lec->deg_line << endl;
	    		deg_lec = deg_lec->next;
	  		}
		}
	}
	output.close();
}
string fileManager::table(string s, int w){	//this function is to write integra.dat files in table format
    stringstream ss, spaces;
    int padding = w - s.size();                 // count excess room to pad
    for(int i=0; i<padding/2; ++i)
        spaces << " ";
    ss << spaces.str() << s << spaces.str();    // format with padding
    if(padding>0 && padding%2!=0)               // if odd #, add 1 space
        ss << " ";
    return ss.str();
}

void fileManager::RepresentationGraphics(string save, int aux, string aprox, int type){
	lines *lector;
	string name;
	if(type==1){
		lector=root;
		name=save+"-VMD.vmd";
	}else{
		lector=root3;
		name=save+"-rho-VMD.vmd";
		save=save+"-rho";
	}

	ofstream output(name.c_str(),ios::trunc);	//output file, name may change.
	string size="1.0",thickness="2", resolution="100", 	color_text="white",color_volume="lime";
	//# graphics delete all
	output << "mol delete all\n"                     ;
	//# graphics top delete all
	output << "draw delete all\n"                    ;
	//# Display settings            
	output << "display projection   Orthographic\n"  ;
	output << "display depthcue     off\n"           ;
	//# axes turn off
	output << "axes location off\n"                  ;
	//# draw basins text
	output << "draw color "<<color_text<<"\n"             ;
	int i=0;
	if ( lector != 0 ) { //Makes sure there is a place to start
		while(i<aux){
		//while(lector->next!=NULL){
			i++;
	  		output <<"graphics top text \""<<lector->c_x<<" "<<lector->c_y<<" "<<lector->c_z<<"\" "<<lector->integral<<" size "<<size<<" thickness "<<thickness<<endl;
	    	lector = lector->next;
	  	}
	}
	//###############################	
	//# isosurface representation
	//# load new molecule 
//	output << "mol new {"<<GetPath()<<"/"<<aprox<<"/"<<save<<".cube" <<"} type cube first 0 last -1 step 1 filebonds 1 autobonds 1 waitfor all\n"  ;
	output << "mol new "<<save<<".cube type cube first 0 last -1 step 1 filebonds 1 autobonds 1 waitfor all\n"  ;
	//# representation of the atoms
	output << "mol delrep 0 top\n"            ;
	output << "mol addrep top\n"              ;
	output << "mol representation CPK 0.500000 0.200000 4000.000000 4000.000000\n"   ;
	output << "mol color Name\n"              ;
	output << "mol selection {all}\n"         ;
	output << "mol material Difusse\n"        ;
	//# add representation of the surface
	output << "mol addrep top\n"              ;
	output << "mol representation Isosurface 0.50000 1 0 0 1 1\n"  ;
	output << "mol color Volume 3\n"          ;
	output << "mol selection {all}\n"         ;
	output << "mol material Transparent\n"    ;
	output << "mol addrep top\n"              ;
	output << "mol selupdate 2 top 0\n"       ;
	output << "mol colupdate 2 top 0\n"       ;
	output << "mol scaleminmax top 2 -0.0000 0.0000\n"  ;
	output << "mol smoothrep top 2 0\n"       ;
	output << "mol drawframes top 2 {now}\n"  ;
	output << "color scale method BGR\n"	     ;
	output.close();
}
void fileManager::WriteMol2(string save, int aux, string aprox){
	//cube= save.cube
	//dat= save+"-integral.txt";
	int size_of_mol2=0, read_line=1, ahora=1;
	string data, xyz="", x, y, z, elementNumber;
	ifstream file;
	file.open((save+".cube").c_str());
	while(file){
		getline(file, data);
		if(read_line>6){
			size_of_mol2++;
			istringstream iss(data);
			ostringstream convert;
			convert << size_of_mol2;
			xyz=xyz+convert.str();
			iss >> elementNumber >> x>>y>>z;
			xyz=xyz+"\t"+x+"\t"+y+"\t"+z+"\n";
			if(ahora==4){
				break;
			}
			ahora++;
		}	
		read_line++;
	}
	file.close();

	cout << xyz<<endl;


	ofstream output((save+".mol2").c_str(),ios::trunc);	//output file, name may change.
	output <<"@<TRIPOS>MOLECULE\n";
	output <<"generated by TAFF\n";
	output <<size_of_mol2+aux<<"\t0\t1\t0\t0\n";
	output <<"SMALL\n";
	output <<"USER_CHARGES\n";
	output <<"****\n";
	output <<"Energy = 0\n";
	output.close();


}
void fileManager::MultiwfnRoute(){	//this function has to be the first one to be call
	ifstream FILE;
	string data;	
	FILE.open("route.txt");
	if (!FILE){
 		cout<<"route.txt not present in current folder, please read manual for more information"<<endl;
 		exit(1);
	}
	while(FILE){
		getline(FILE,data);
		size_t f1=data.find("Multiwfn = ");

		if(f1!=string::npos && data.substr(0,1)!="#"){//the data line has "Mulltiwfn = " and doesn't start with #
			istringstream iss(data.substr(data.find("=")+1));
			iss >> Multiwfn_Route;
		}
	}
	FILE.close();
	//cout <<Multiwfn_Route<<endl; //this shows the Multiwfn program route.
}
string fileManager::currentDateTime() {	//again to write time and date in -integral.dat files
    time_t     now = time(0);
    struct tm  tstruct;
    char       buf[80];
    tstruct = *localtime(&now);
    strftime(buf, sizeof(buf), "%Y-%m-%d.%X", &tstruct);

    return buf;
}

void fileManager::PrintPcenOW(long double Occ[], int OccSize, long double Unocc[] , int UnoccSize, long double sumFMinus, long double sumFPlus, string name){
//"f+orbitalweighted-"+name
	//f-
	string save="f-orbitalweighted-"+name+"-integral.txt";
	ofstream output(save.c_str(),ios::app);	//output file, name may change.

	//contribucion 
	output <<"\n·······································································································\n";
	output <<"ORBITAL REACTIVITY CONTRIBUTION\n·······································································································\nOrbital\tWeight\tContribution (%)"<<endl;
	//add for 
	output << "HOMO\t"<< (Occ[OccSize-1]) << "\t"<<((Occ[OccSize-1])*100)/sumFMinus<<"%"<<endl;
	for (int i = 1; i < OccSize; ++i){
		output << "HOMO-"<<i<<"\t"<< (Occ[OccSize-i-1]) << "\t"<<((Occ[OccSize-i-1])*100)/sumFMinus<<endl;
	}
	//f+
	/*
	save="f+orbitalweighted-"+name+"-integral.txt";
	ofstream output2(save.c_str(),ios::app);	//output file, name may change.
	output2 << "\nLUMO\t"<< (Unocc[0]) << "\t"<<(Unocc[0]*100)/sumFPlus<<"%"<<endl;
	for (int i = 1; i <= UnoccSize; ++i){
		output2 << "LUMO+"<<i<<"\t"<< (Unocc[i]) << "\t"<<(Unocc[i]*100)/sumFPlus<<"%"<<endl;
	}
	*/
	//cout << "sumFPlus "<<sumFPlus<< " sumFMinus "<<sumFMinus <<endl;
	//profit
}
