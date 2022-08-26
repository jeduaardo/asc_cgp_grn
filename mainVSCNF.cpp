#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <time.h>
#include "cgp.h"
#include <string.h>
#include <map>
#include<cstdlib>
#include <sstream> 
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>


using namespace std;


void gera_rede_saida(int ni, int nTests, vector<vector<int>> OTIMIZADOS_LOCAL, char *dir, char *pt){


	ofstream rankedEdges;
  


  string directory = dir;

  string separa = "/";

	string filename = "rankedEdges_";

  string pseudotime = pt;

  string extensao = ".csv";

  string final = directory + separa + filename + pt + extensao;

	rankedEdges.open(final);
	map<int, string> dicionario;


  //VSC Genes
  dicionario[0] = "Dbx1";
  dicionario[1] = "Dbx2";
  dicionario[2] = "Irx3";
  dicionario[3] = "Nkx22";
  dicionario[4] = "Nkx61";
  dicionario[5] = "Nkx62";
  dicionario[6] = "Olig2";
  dicionario[7] = "Pax6";
  dicionario[8] = "Dbx1";
  dicionario[9] = "Dbx2";
  dicionario[10] = "Irx3";
  dicionario[11] = "Nkx22";
  dicionario[12] = "Nkx61";
  dicionario[13] = "Nkx62";
  dicionario[14] = "Olig2";
  dicionario[15] = "Pax6";



	int metade = ni/2;



  for(int i = 0; i < OTIMIZADOS_LOCAL.size(); i++){
    vector<int> vetor_local = {};
    for(int j = 0; j < ni; j++){
      int mycount = count(OTIMIZADOS_LOCAL[i].begin(), OTIMIZADOS_LOCAL[i].end(), j);
      vetor_local.push_back(mycount);
    }
    for(int k = 0; k < vetor_local.size(); k++){
      if(vetor_local[k] != 0){
        rankedEdges << dicionario[k] << "\t" << dicionario[i] << "\t";
        float conta = (float)vetor_local[k] / (float)nTests;
        if(k >= metade){
          rankedEdges << "-" << conta;
        }
        else{
          rankedEdges << conta;
        }
        rankedEdges << "\n";
      }
    }
  }

}



void printFirstSolution(Individual * p, string Letra, string test, char *dir, char *pt){
	ofstream arquivo_factivel;
    string directory = dir;
    string pseudotime = pt;
    string filename = directory + "/" + Letra + "_" + test + "_" + pseudotime + "_feasible.txt";
    arquivo_factivel.open(filename);
    arquivo_factivel << "Test Number: " << test << endl;
    for(int outputs = 0; outputs < p->no; outputs++){
    	arquivo_factivel << "Output " << outputs << endl;
    	sort(p->phenotype[outputs].begin(), p->phenotype[outputs].end());
    	vector<int> phenotype_local = {};

    	for(int nodes = 0; nodes < p->phenotype[outputs].size(); nodes++){

			auto search = find(phenotype_local.begin(), phenotype_local.end(), p->phenotype[outputs][nodes]);
			if(search == phenotype_local.end()){
				phenotype_local.push_back(p->phenotype[outputs][nodes]);
			}	


    	}


    	for(int actives = 0; actives < phenotype_local.size(); actives++){
        	int currentActive = phenotype_local[actives];
        	vector<int> gene = p->genotype[currentActive];
        	arquivo_factivel << "NODE: " << currentActive << " - ";
        	for(int i = 0; i < gene.size(); i++){
          		arquivo_factivel << gene[i] << "  ";
        	}
        	arquivo_factivel << "\n";
    	}
    }
    arquivo_factivel.close();
}



void printFinalSolution(vector<vector<int>> genotype, string Letra, string test, int ni_c, int no_c, int nc_c, int nr_c, int lb_c, char *dir, char *pt){
	Individual * local = new Individual(ni_c, no_c, nc_c, nr_c, lb_c, genotype);
  	for(int j = 0; j < local->no; j++){
    	int currentOutput = local->size + local->ni + j;
      	local->phenotype[j].push_back(currentOutput);
      	local->phenotype[j].push_back(local->genotype[currentOutput][0]);
      	local->getActiveNodes(j, local->genotype[currentOutput][0]);
  	}
    string directory = dir;
    string pseudotime = pt;
  	ofstream arquivo_final;
  	string filename = directory + "/" + Letra + "_" + test + "_" + pt + "_optimized.txt";
  	arquivo_final.open(filename);
  	for(int outputs = 0; outputs < local->no; outputs++){
    	arquivo_final << "Output " << outputs << endl;
      	sort(local->phenotype[outputs].begin(), local->phenotype[outputs].end());
      
    	vector<int> phenotype_local = {};

    	for(int nodes = 0; nodes < local->phenotype[outputs].size(); nodes++){

			auto search = find(phenotype_local.begin(), phenotype_local.end(), local->phenotype[outputs][nodes]);
			if(search == phenotype_local.end()){
				phenotype_local.push_back(local->phenotype[outputs][nodes]);
			}	


    	}

      	for(int actives = 0; actives < phenotype_local.size(); actives++){
        	int currentActive = phenotype_local[actives];
          	vector<int> gene = local->genotype[currentActive];
          	arquivo_final << "NODE: " << currentActive << " - ";
          	for(int i = 0; i < gene.size(); i++){
              	arquivo_final << gene[i] << "  ";
          	}
      		arquivo_final << "\n";
      	}
  	}
  	local->countLE3();
  	arquivo_final.close();
}



void printNonFeasibleSolution(vector<vector<int>> genotype, string Letra, string test, int ni_c, int no_c, int nc_c, int nr_c, int lb_c, char *dir, char *pt){
	Individual * local = new Individual(ni_c, no_c, nc_c, nr_c, lb_c, genotype);
  	for(int j = 0; j < local->no; j++){
    	int currentOutput = local->size + local->ni + j;
      	local->phenotype[j].push_back(currentOutput);
      	local->phenotype[j].push_back(local->genotype[currentOutput][0]);
      	local->getActiveNodes(j, local->genotype[currentOutput][0]);
  	}

    string directory = dir;
    string pseudotime = pt;
    ofstream arquivo_final;
    string filename = directory + "/" + Letra + "_" + test + "_" + pt + "_nonFeasible.txt";

  	arquivo_final.open(filename);
  	for(int outputs = 0; outputs < local->no; outputs++){
    	arquivo_final << "Output " << outputs << endl;
      	sort(local->phenotype[outputs].begin(), local->phenotype[outputs].end());
      
    	vector<int> phenotype_local = {};

    	for(int nodes = 0; nodes < local->phenotype[outputs].size(); nodes++){

			auto search = find(phenotype_local.begin(), phenotype_local.end(), local->phenotype[outputs][nodes]);
			if(search == phenotype_local.end()){
				phenotype_local.push_back(local->phenotype[outputs][nodes]);
			}	


    	}

      	for(int actives = 0; actives < phenotype_local.size(); actives++){
        	int currentActive = phenotype_local[actives];
          	vector<int> gene = local->genotype[currentActive];
          	arquivo_final << "NODE: " << currentActive << " - ";
          	for(int i = 0; i < gene.size(); i++){
              	arquivo_final << gene[i] << "  ";
          	}
      		arquivo_final << "\n";
      	}
  	}
  	arquivo_final.close();
}


int main(int argc, char *argv[ ])
{ 


      char discretizationMethod[50];
      char problemConfig[50];
      char pseudotime[5];
      char tabelaverdade[70] = "inputs/";
      char texto[10] = "";




      strcpy(discretizationMethod, argv[1]);
      strcpy(problemConfig, argv[2]);
      strcpy(pseudotime, argv[3]);
      strcat(tabelaverdade, argv[4]);
      strcat(tabelaverdade, texto);
     
      cout << tabelaverdade << endl;

      char separa[2] = "/";

      char DM[50];
      strcpy(DM, discretizationMethod);



      strcat(DM, separa);
      strcat(DM, problemConfig);



      mkdir(discretizationMethod, S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
      mkdir(DM, S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);


   vector<vector<int>> dAllOutputs = {};
   string myText;
   ifstream MyReadFile(tabelaverdade);
   int linha = 0;
   while (getline (MyReadFile, myText)) {
      dAllOutputs.push_back({});
      for(int i = 0; i < myText.size(); i++){
         int valor_atual = myText[i];
         int ivalor_atual = valor_atual - '0';
         dAllOutputs[linha].push_back(ivalor_atual);
      }

      //cout << myText;
      linha += 1;
   }
   MyReadFile.close();
  

   //cout << dAllOutputs[0].size() << endl;
   map<string, int> intparam;
   map<string, string> strparam;
   map<int, string> alfabeto;
   intparam["nTests"] = 5;
   intparam["nEvaluations"] = 50000;
   intparam["initialSeed"] = 16;
   intparam["mu"] = 1;
   intparam["lambda"] = 4;
   intparam["ni"] = 16;
   intparam["no"] = 1;
   intparam["nc"] = 100;
   intparam["nr"] = 1;
   intparam["lb"] = 100;
   alfabeto[0] = "A"; 
   alfabeto[1] = "B";
   alfabeto[2] = "C";
   alfabeto[3] = "D";
   alfabeto[4] = "E";
   alfabeto[5] = "F";
   alfabeto[6] = "G";
   alfabeto[7] = "H";
   alfabeto[8] = "I";
   alfabeto[9] = "J";
   alfabeto[10] = "K";
   alfabeto[11] = "L";
   alfabeto[12] = "M";
   alfabeto[13] = "N";
   alfabeto[14] = "O";
   alfabeto[15] = "P";
   alfabeto[16] = "Q";
   alfabeto[17] = "R";
   alfabeto[18] = "S";
   alfabeto[19] = "T";
   alfabeto[20] = "U";
   alfabeto[21] = "V";
   alfabeto[22] = "W";
   alfabeto[23] = "X";
   alfabeto[24] = "Y";
   alfabeto[25] = "Z";                        

   vector<vector<int>> OTIMIZADOS = {};

   for(int saida = 0; saida < dAllOutputs.size(); saida++){
   	OTIMIZADOS.push_back({});
   		vector<vector<int>> dOutput = {dAllOutputs[saida]};
   		for(int test = 0; test < intparam["nTests"]; test++){
   			srand(intparam["initialSeed"]+test);
   			//cout << "Initializing Test Number " << test << "/" << intparam["nTests"]-1 << " for variable " << alfabeto[saida] <<endl;
   			
   			Population * k = new Population(intparam["mu"],intparam["lambda"], dOutput); 
   			k->clearPopulation();

   			for (int i=0; i<intparam["mu"]+intparam["lambda"]; i++){
      			k->Pop.push_back(new Individual(intparam["ni"], intparam["no"], intparam["nc"], intparam["nr"], intparam["lb"], {}));
      			k->Pop[i]->code_genotype();
      			for(int j = 0; j < k->Pop[i]->no; j++){
         			int currentOutput = k->Pop[i]->size + k->Pop[i]->ni + j;
    	 			k->Pop[i]->phenotype[j].push_back(currentOutput);
    	 			k->Pop[i]->phenotype[j].push_back(k->Pop[i]->genotype[currentOutput][0]);
    	 			k->Pop[i]->getActiveNodes(j, k->Pop[i]->genotype[currentOutput][0]);
      			}
	   
	  			k->Pop[i]->evalIndividual(k->desiredOutputs);
		    }


    		vector<int> all_fitness;
 
    		for(int i =0; i<intparam["mu"]+intparam["lambda"];i++){
      			all_fitness.push_back(k->Pop[i]->fitness);
      			if(k->Pop[i]->fitness == (pow(2, intparam["ni"]/2)) * k->Pop[i]->no){
        			//cout << "Feasible Solution Found on Initial Population" << endl;
        			string nome2;
        			if (test == 0){
          				nome2 = "0";
        			}
        			if (test == 1){
          				nome2 = "1";
        			}
        			if (test == 2){
          				nome2 = "2";
        			}
        			if (test == 3){
          				nome2 = "3";
        			}
        			if (test == 4){
          				nome2 = "4";
        			}                                
        			printFirstSolution(k->Pop[i], alfabeto[saida], nome2, DM, pseudotime);
        			k->parentUpdateOptimize1(k->Pop[i]->genotype, k->Pop[i]->fitness, k->Pop[i]->LE);
        			k->feasible = 1;        
      			}
    		}

    		k->parentUpdateIndividual(k->Pop[max_element( all_fitness.begin(), all_fitness.end()) - all_fitness.begin()]->genotype, k->Pop[max_element( all_fitness.begin(), all_fitness.end()) - all_fitness.begin()]->fitness, intparam["ni"], intparam["no"], intparam["nc"], intparam["nr"], intparam["lb"]); // ATUALIZA CORRETAMENTE O BEST

    		int evaluations = 0;
    		while(evaluations < intparam["nEvaluations"]){
    			k->clearPopulation();

   				for (int i=0; i<intparam["lambda"]; i++){
        			k->Pop.push_back(new Individual(intparam["ni"], intparam["no"], intparam["nc"], intparam["nr"], intparam["lb"], k->best));
        
        			for(int j = 0; j < k->Pop[i]->no; j++){
            			int currentOutput = k->Pop[i]->size + k->Pop[i]->ni + j;
            			k->Pop[i]->phenotype[j].push_back(currentOutput);
            			k->Pop[i]->phenotype[j].push_back(k->Pop[i]->genotype[currentOutput][0]);
            			k->Pop[i]->getActiveNodes(j, k->Pop[i]->genotype[currentOutput][0]);
         			}
              k->Pop[i]->mutateSAM();  

        			
         			for(int j = 0; j < k->Pop[i]->no; j++){
            			int currentOutput = k->Pop[i]->size + k->Pop[i]->ni + j;
            			k->Pop[i]->phenotype[j].push_back(currentOutput);
            			k->Pop[i]->phenotype[j].push_back(k->Pop[i]->genotype[currentOutput][0]);
            			k->Pop[i]->getActiveNodes(j, k->Pop[i]->genotype[currentOutput][0]);
         			}       
      
        			k->Pop[i]->evalIndividual(k->desiredOutputs);
        			evaluations += 1;
        			if(k->feasible == 1){
          				k->Pop[i]->countLE();
        			}
    			}


    			all_fitness.clear();
    			int max_fitness = pow(2, intparam["ni"]/2) * k->Pop[0]->no;


    			for(int i =0; i<intparam["lambda"];i++){
      				all_fitness.push_back(k->Pop[i]->fitness);
      
      
      				if(k->Pop[i]->fitness == max_fitness && k->feasible == 0){
        				//cout << "Feasible Solution Found with " << evaluations << " evaluations." << endl;
        				string nome;
        				if (test == 0){
          					nome = "0";
        				}
        				if (test == 1){
          					nome = "1";
        				}
        				if (test == 2){
          					nome = "2";
        				}
        				if (test == 3){
          					nome = "3";
        				}
        				if (test == 4){
          					nome = "4";
        				}                                
        				printFirstSolution(k->Pop[i], alfabeto[saida], nome, DM, pseudotime);
        				k->parentUpdateOptimize1(k->Pop[i]->genotype, k->Pop[i]->fitness, k->Pop[i]->LE);
        				k->feasible = 1;
        				break;
      				}
    			}

    			int greaterFitness = k->Pop[max_element( all_fitness.begin(), all_fitness.end()) - all_fitness.begin()]->fitness;

    			if(k->feasible == 0){
      				if(greaterFitness >= k->currentFitness){
        				k->parentUpdateIndividual(k->Pop[max_element( all_fitness.begin(), all_fitness.end()) - all_fitness.begin()]->genotype, k->Pop[max_element( all_fitness.begin(), all_fitness.end()) - all_fitness.begin()]->fitness, intparam["ni"], intparam["no"], intparam["nc"], intparam["nr"], intparam["lb"]); // ATUALIZA CORRETAMENTE O BEST
      				}

			    }
    			else{

      				k->parentUpdateOptimize();

    			}    

    			if(evaluations % 1000 == 0){
      				cout << "Evaluations: " << evaluations << " Fitness: " << k->currentFitness << endl;
      				cout << "Logic Elements " << k->currentLE << endl;
            int tg = 0;
    			}

    		} // Fim Processo Geracional WHILE
    		cout << "End of Test Number " << test << "/" << intparam["nTests"]-1 << " for variable " << alfabeto[saida] <<endl;
    		string nome2;
        	if (test == 0){
          		nome2 = "0";
        	}
        	if (test == 1){
          		nome2 = "1";
        	}
        	if (test == 2){
          		nome2 = "2";
        	}
        	if (test == 3){
          		nome2 = "3";
        	}
        	if (test == 4){
          		nome2 = "4";
        	} 
        	if(k->feasible == 1){
	    		printFinalSolution(k->best, alfabeto[saida], nome2, intparam["ni"], intparam["no"], intparam["nc"], intparam["nr"], intparam["lb"], DM, pseudotime);


				Individual * local = new Individual(intparam["ni"], intparam["no"], intparam["nc"], intparam["nr"], intparam["lb"], k->best);
  				for(int j = 0; j < local->no; j++){
    				int currentOutput = local->size + local->ni + j;
      				local->phenotype[j].push_back(currentOutput);
      				local->phenotype[j].push_back(local->genotype[currentOutput][0]);
      				local->getActiveNodes(j, local->genotype[currentOutput][0]);
  				}
  	
  
  				for(int outputs = 0; outputs < local->no; outputs++){
      				sort(local->phenotype[outputs].begin(), local->phenotype[outputs].end());
      
    				vector<int> phenotype_local = {};

    				for(int nodes = 0; nodes < local->phenotype[outputs].size(); nodes++){

						auto search = find(phenotype_local.begin(), phenotype_local.end(), local->phenotype[outputs][nodes]);
						if(search == phenotype_local.end()){
							phenotype_local.push_back(local->phenotype[outputs][nodes]);
						}	
    				}

      				for(int actives = 0; actives < phenotype_local.size(); actives++){
        				int currentActive = phenotype_local[actives];
        				if (currentActive < intparam["ni"]){
        					OTIMIZADOS[saida].push_back(currentActive);
        				}
      				}
  				}


    		}
    		else{
    			printNonFeasibleSolution(k->best, alfabeto[saida], nome2, intparam["ni"], intparam["no"], intparam["nc"], intparam["nr"], intparam["lb"], DM, pseudotime);
    		}
    		k->setStandard();
		} // Fim nTests
} //Fim nSaidas
gera_rede_saida(intparam["ni"], intparam["nTests"], OTIMIZADOS, DM, pseudotime);
return 0;
} // Fim main
