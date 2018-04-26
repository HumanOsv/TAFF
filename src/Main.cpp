#include <iostream>
#include <vector>
#include "menu.h"
#include <time.h>
#include <sys/time.h>
std::string currentDateTime() {	//simple function to get actual date and time
    time_t     now = time(0);
    struct tm  tstruct;
    char       buf[80];
    tstruct = *localtime(&now);
    strftime(buf, sizeof(buf), "%Y-%m-%d.%X", &tstruct);
    return buf;
}//Authors: O. Yañez, D.Inostroza,  R. Pino-Rios   ,  A.  Vásquez-Espinal,  P. Fuentealba    & W. Tiznado

int main(int argc, char* argv[]) { 
	//Main function, it send arguments information to menu class
	time_t before=time(0);
	std::cout << "\n" ;
	std::cout << "\t  __________    ____    ________   ________\n" ;
	std::cout << "\t  \\ __   ___\\  / __ \\  |  _____/  |  _____/\n" ;
	std::cout << "\t      | |     / /  \\ \\ | |        | |\n"       ;	
	std::cout << "\t      | |     | |__| | | |___     | |___\n"    ;	
	std::cout << "\t      | |     |  __  | |  __/     |  __/\n"    ;	
	std::cout << "\t      | |     | |  | | | |        | |\n"       ;	
	std::cout << "\t      |_|     |_|  |_| |_|        |_|\n"       ;	
	std::cout <<  "" ;	
	std::cout << "\t# # # # # # # # # # # # # # # # # # # # # # # # # #\n" ;	
	std::cout << "\t# Topological Analysis of Fukui Functions (TAFF)  #\n" ;	  
	std::cout << "\t# DATE: "<<currentDateTime()<<"\t\t          #\n";	
	std::cout << "\t# TiznadoLab\t\t\t\t\t  #\n";
	std::cout << "\t# # # # # # # # # # # # # # # # # # # # # # # # # #\n" ;	 
	std::cout <<  "\n" ;	
	std::cout <<  "\n" ;	

	std::vector<menu> men;
	menu *chain=new menu();

	chain->ArgumentsValidation(argc, argv);
	chain->Mainmenu();
	delete chain;
	time_t after= time(0);
	int sad=after-before;
	std::cout << "\nProgram finished correctly"<<std::endl;
	std::cout << "CPU TIME: "<<sad <<" seconds"<<std::endl;
	

	return 0;


}