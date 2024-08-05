/*
* NAME	: LUM FU YUAN
* MATRIC: B032110251
* GROUP	: BITI S1G2
*/

#include <iostream>
#include <iomanip>
#include <cstdlib>
#include <ctime>
#include <string>
#include <fstream>
#include <algorithm>
#include <vector>

using namespace std;

//##################################  constants  #########################################
const int MAP_SIZE = 10;
const int GENE = MAP_SIZE;
const int POP_SIZE = 200;
const int BUDGET = 1500;
//const int PENALTY = 1500; 
const int MAP_AREA = MAP_SIZE * MAP_SIZE;
const double CCTV_PRICE[] = {/*A*/75,/*B*/200};

const double CROSSOVER_PROBABILITY = 0.8;
const double MUTATION_PROBABILITY = 0.2;
const int MAXIMUM_GENERATION = 55;

// #################################  data storage  ####################################

string chromosome[POP_SIZE][GENE];
string parents[2][GENE];
string children[2][GENE];
string NewChromosome[POP_SIZE][GENE];

int mapValue[MAP_SIZE][MAP_SIZE];
double fitness_value[POP_SIZE];		
int countNewChromosome = 0;

string bestChromosome[GENE];
float bestFitness = -1;
double avgFitness;
// #################################  file storage  ####################################
ofstream avgFitnessFile;
ofstream BestFitnessFile;
ofstream BestChromosomeFile;
ofstream BestSolutionFile; 

//#####################################  function prototype  ##############################
void initialMap();
void initializePopulation();
void printChromosome();
void evaluateChromosome();
void parentSelection();
void crossover();
void mutation();
void survivalSelection();
void copyChromosome();
void calculateAverageFitness();
void recordBestFitness();
//for best solution
void GenerateBestSolutionFile(int Gen);


void displayMapValue();

int main() {
	//open file
	avgFitnessFile.open("AvgFinessFile.txt");
	BestFitnessFile.open("BestFitnessFile.txt");
	BestChromosomeFile.open("BestChromosomeFile.txt");
	BestSolutionFile.open("BestSolution.csv");

	BestSolutionFile << "Generation,Average Fitness,Best Fitness,,,Best Chromosome" << endl; // for Excel header

	cout << "Initialize population:\n";
	initializePopulation();
	printChromosome();

	//loop generation
	for (int G = 0; G < MAXIMUM_GENERATION; G++) {
		countNewChromosome = 0;
		cout << "\n========================== Generation " << G + 1 << "===========================\n";

		cout << "\nEvaluate chromosome\n";
		evaluateChromosome();
		calculateAverageFitness();
		recordBestFitness();
		GenerateBestSolutionFile(G);


		for (int i = 0; i < POP_SIZE / 2; i++) {
			cout << "\nParent Selection\n";
			parentSelection();

			cout << "\nCrossover Method\n";
			crossover();

			cout << "\nMutation Method\n";
			mutation();

			cout << "\n\nSurvival Selection Method\n";
			survivalSelection();
		}
		cout << "\n\nCopy Chromosome....\n";
		copyChromosome();
	}

	avgFitnessFile.close();
	BestFitnessFile.close();
	BestChromosomeFile.close();
	BestSolutionFile.close();
	return 0;
}

void initialMap() {
	// intital map value=0;
	for (int r = 0; r < MAP_SIZE; r++) {
		for (int c = 0; c < MAP_SIZE; c++) {
			mapValue[r][c] = 0;
		}
	}
}

void initializePopulation() {
	srand(time(NULL));
	int randNum;
	int allele1;
	string allele2;
	const string TYPE_CCTV[] = { "A", "B" };

	for (int c = 0; c < POP_SIZE; c++) {
		for (int g = 0; g < GENE; g++) {
			allele1 = rand() % (MAP_SIZE+1);      //remainder must <MAP_SIZE+1 
			if (allele1 != 0) {
				randNum = rand() % 2;	
				allele2 = TYPE_CCTV[randNum];	//random select type CCTV(A,B)
			}
			else {
				allele2 = "0";
			}
			chromosome[c][g] = to_string(allele1) + allele2;
		}
	}
}

void printChromosome() {
	for (int c = 0; c < POP_SIZE; c++) {
		cout << "Chromosome " << c + 1 << "\t: ";
		for (int g = 0; g < GENE; g++) {
			cout << setw(3) << chromosome[c][g] << " ";
		}
		cout << endl;
	}
}

void displayMapValue(int type) {
	if (type == 1) { //for console
		for (int c = 0; c < MAP_SIZE; c++) {
			cout << "row        " << c + 1 << "\t: ";
			for (int g = 0; g < MAP_SIZE; g++) {
				cout << setw(3) << mapValue[c][g] << " ";
			}
			cout << endl;
		}
	}
	else if (type == 2) {// for Excel
		BestSolutionFile << endl;
		for (int c = 0; c < MAP_SIZE; c++) {
			BestSolutionFile << ",,,,row" << c + 1 << ":,";
			for (int g = 0; g < MAP_SIZE; g++) {
				BestSolutionFile << mapValue[c][g] << ",";
			}
			BestSolutionFile << endl ;
		}
	}
}

void evaluateChromosome() {
	int NumCCTV;
	double SumVisionValue,SumCCTVPrice;
	double f1, f2;
	int allele1;
	string allele2;

	for (int c = 0; c < POP_SIZE; c++) {
		NumCCTV = 0;
		f1 = f2 = SumVisionValue = SumCCTVPrice = 0.0;
		initialMap();
		
		for (int g = 0; g < GENE; g++) {

			if (chromosome[c][g] != "00") {
				NumCCTV++;

				//sample = "3B" , allele1=3, allele2=B 
				if (chromosome[c][g].length() < 3) {
					allele1 = stoi(chromosome[c][g].substr(0, 1));	//Starting position = 0, length = 1 (get num after the position)
					allele2 = chromosome[c][g].substr(1, 1);		//Starting position = 1, length = 1 (get num after the position)
				}
				else {	// need handle when 2 decimal sample = "10A" , allele1=10, allele2=A
					allele1 = stoi(chromosome[c][g].substr(0, 2));	//Starting position = 0, length = 1 (get num after the position)
					allele2 = chromosome[c][g].substr(2, 1);		//Starting position = 1, length = 1 (get num after the position)
				}

				int y = allele1; 
				int x = g+1;	//(1-10)

				//set map vision area (x=c, y=r)
				if (allele2 == "A") {
					//A CCTV
					int r = y - 2;
					int c = x - 2;
					for (int i = r; i < r + 3; i++) {
						for (int j = c; j < c + 3; j++) {
							if ((i >= 0 && i < MAP_SIZE) && (j >= 0 && j < MAP_SIZE)) {
								mapValue[i][j] =1;
							}
						}
					}
				}
				else if (allele2 == "B") {
					//B CCTV
					int r = y - 3;
					int c = x - 3;
					for (int i = r; i < r + 5; i++) {
						for (int j = c; j < c + 5; j++) {
							if ((i >= 0 && i < MAP_SIZE) && (j >= 0 && j < MAP_SIZE)) {
								mapValue[i][j] = 1;
							}
						}
					}
				}
				else {
					cout << "something went wrong at evaluateChromosome";
					exit(0);
				}
				
				// calculate total CCTV price

				if (allele2 == "A") {
					SumCCTVPrice += CCTV_PRICE[0];
				}
				else if (allele2 == "B") {
					SumCCTVPrice += CCTV_PRICE[1];
				}
				else {
					cout << "something went wrong at calculate total CCTV price";
					exit(0);
				}

			}//end if(chromosome[c][g] != "00")

		}//end one chromosome loop

		//display each chromosome vision area 2D
		
		//for console to display the map of Coverage of the CCTV 
		/*
				//for console
				cout << " x:"<<x << " y:" << y<<" type"<< allele2 << endl;
				displayMapValue(1);
				cout << endl;
		*/
		
		


		// calculate total vision area
		for (int r = 0; r < MAP_SIZE; r++) {
			for (int c = 0; c < MAP_SIZE; c++) {
				SumVisionValue += mapValue[r][c];
			}
		}

		// calculate fitness 1
		f1 = SumVisionValue / (MAP_SIZE * MAP_SIZE);

		// calculate fitness 2			
		f2 = (max(static_cast<double>(BUDGET), SumCCTVPrice) - SumCCTVPrice) / max(static_cast<double>(BUDGET), SumCCTVPrice);


		// fitness value
		fitness_value[c] = (f1 * 0.5) + (f2 * 0.5);
		cout << "chromosome: " << c + 1 << "\tNumCCTV: "<< NumCCTV << "\tSumVisionValue: " << SumVisionValue << "\tSumCCTVPrice: " << SumCCTVPrice << "\tFitness1: " << f1 << "\tFitness2: " << f2 << "\tTotalFitness: " << fitness_value[c] << endl;
		
	}//end one population loop

}

void parentSelection() {
	int player1, player2;
	int indexParents[2];

	// tournament selection

	do {
		//1. For both parents
		for (int p = 0; p < 2; p++) {
			//pick a random number to be the index for player 1
			player1 = rand() % POP_SIZE;
			do {
				//pick another random number to be the index for player 2
				player2 = rand() % POP_SIZE;
			} while (player1 == player2); //avoid player 1 same with player 2

			cout << "Player " << player1 + 1 << " vs Player " << player2 + 1 << endl;

			// determine indexParents
			if (fitness_value[player1] >= fitness_value[player2]) {
				indexParents[p] = player1;
			}
			else {
				indexParents[p] = player2;
			}

			cout << "\n\t Player's index: " << player1+1 << " vs " << player2+1;
			cout << "\n\t Fitness: " << fitness_value[player1] << " vs " << fitness_value[player2];
			cout << "\n\t Winner " << p + 1 << " : " << indexParents[p]+1 << endl << endl;
		} // end of tournament

	} while (indexParents[0] == indexParents[1]); // restart tournament if same indexParents

	//Print Parent
	for (int p = 0; p < 2; p++) {

		cout << "\tParent " << p + 1 << "\t";

		for (int i = 0; i < GENE; i++) {
			parents[p][i] = chromosome[indexParents[p]][i];
			cout << setw(3) << parents[p][i] << " ";
		}
		cout << "\n";
	}
}

void crossover() {
	double randNum;
	//=========================================================== uniform crossover ====================================
	// Copy both parent’s chromosome to children chromosomes
	for (int c = 0; c < 2; c++) {
		for (int g = 0; g < GENE; g++) {
			children[c][g] = parents[c][g];
		}
	}

	// Generate a random number from 0-1. Make sure it is a real value data type
	randNum = (rand() % 11) / 10.0;

	cout << "\n\tProbablility : " << CROSSOVER_PROBABILITY;

	// If randNum is less than crossover probability
	if (randNum <= CROSSOVER_PROBABILITY) {
		cout << "\n\tCrossover happens";

		vector<int> crossoverIndex;
		// Perform uniform crossover for each gene
		for (int g = 0; g < GENE; g++) {
			if (rand() % 2 == 0) {
				crossoverIndex.push_back(g);
				swap(children[0][g], children[1][g]);
			}
		}
		cout << "\n\tCrossover points: ";
		for (int point : crossoverIndex) {
			cout << point << " ";
		}

	}
	else {
		cout << "\n\tCrossover does not happen";
	}

	cout << "\n";
	cout << "\n\tindex:            0   1   2   3   4   5   6   7   8   9";

	// Print children 1 & 2
	for (int c = 0; c < 2; c++) {
		cout << "\n";
		cout << "\tchildren " << c + 1 << " : " << "\t";
		for (int g = 0; g < GENE; g++) {
			cout << setw(3) << children[c][g] << " ";
		}
	}
	cout << "\n";

	/*
	double randNum;
	int n = 5; // Number of crossover points

	//=========================================================== n-point crossover =====================================
	// Copy both parent’s chromosome to children chromosomes
	for (int c = 0; c < 2; c++) {
		for (int g = 0; g < GENE; g++) {
			children[c][g] = parents[c][g];
		}
	}

	// Generate a random number from 0-1. Make sure it is a real value data type
	randNum = (rand() % 11) / 10.0;

	cout << "\n\tProbablility : " << CROSSOVER_PROBABILITY;

	// If randNum is less than crossover probability
	if (randNum <= CROSSOVER_PROBABILITY) {
		cout << "\n\tCrossover happens";

		// Generate n random crossover points
		vector<int> crossoverPoints;
		for (int i = 0; i < n; i++) {
			int point = rand() % GENE;
			crossoverPoints.push_back(point);
		}

		// Sort the crossover points in ascending order
		sort(crossoverPoints.begin(), crossoverPoints.end());

		cout << "\n\tCrossover points: ";
		for (int point : crossoverPoints) {
			cout << point << " ";
		}

		// Perform crossover at the selected points
		int parentIndex = 0;
		int crossoverIndex = 0;
		for (int g = 0; g < GENE; g++) {
			if (g == crossoverPoints[crossoverIndex]) {
				// Swap parent index when a crossover point is reached
				parentIndex = (parentIndex + 1) % 2;
				crossoverIndex++;

				if (crossoverIndex >= n) {
					// Break the loop if all crossover points are handled
					break;
				}
			}

			// Assign genes from the selected parent to the children
			children[0][g] = parents[parentIndex][g];
			children[1][g] = parents[(parentIndex + 1) % 2][g];
		}
	}
	else {
		cout << "\n\tCrossover does not happen";
	}

	cout << "\n";
	cout << "\n\tindex:            0   1   2   3   4   5   6   7   8   9";

	// Print children 1 & 2
	for (int c = 0; c < 2; c++) {
		cout << "\n";
		cout << "\tchildren " << c + 1 << " : " << "\t";
		for (int g = 0; g < GENE; g++) {
			cout << setw(3) << children[c][g] << " ";
		}

	}
	cout << "\n";
	*/
	/*
	double randNum;
	//=========================================================== one point crossover  ==================================
	// Copy both parent’s chromosome to children chromosomes
	for (int c = 0; c < 2; c++) {
		for (int g = 0; g < GENE; g++) {
			children[c][g] = parents[c][g];
		}
	}

	//Generate a random number from 0-1. Make sure it is real value data type
	randNum = (rand() % 11) / 10.0;

	cout << "\n\tProbablility : " << CROSSOVER_PROBABILITY;

	//If randNum less than crossover probability
	if (randNum <= CROSSOVER_PROBABILITY) {
		cout << "\n\tCrossover happen";

		//generate a random crossover point
		int crossoverPoint = rand() % GENE;

		//Crossover parent 1 and parent 2 to produce the children
		cout << "\n\tCrossover point: " << crossoverPoint;

		for (int g = crossoverPoint; g < GENE; g++) {
			children[0][g] = parents[1][g];
			children[1][g] = parents[0][g];
			//			swap(children[0][g], children[1][g]);
		}
	}
	else {
		cout << "\n\tCrossover no happen";
	}

	cout << "\n";
	cout << "\n\tindex:            0   1   2   3   4   5   6   7   8   9";

	//Print children 1 & 2
	for (int c = 0; c < 2; c++) {
		cout << "\n";
		cout << "\tchildren " << c + 1 << " : " << "\t";
		for (int g = 0; g < GENE; g++) {
			cout << setw(3) << children[c][g] << " ";
		}

	}
	cout << "\n";
	*/
	
}

void mutation() {
	//=========================================================== swap mutaion  ==================================
	cout << "\n\tProbablility : " << MUTATION_PROBABILITY<<endl;
	float randNum, mutationBit;

	//For both children
	for (int c = 0; c < 2; c++) {

		//Generate number from 0-1 (real values)
		randNum = (rand() % 11) / 10.0;

		//If randNum less than mutation probability
		if (randNum <= MUTATION_PROBABILITY) {
			cout << "\n\tMutation happen in children " << c + 1;

			//generate a mutation bit 1
			int mutaionBit1 = rand() % GENE;
			//generate a mutation bit 2
			int mutaionBit2 = rand() % GENE;

			//Print the mutation bit
			cout << "\n\tThe first mutation bit:  " << setw(2) << mutaionBit1 << endl;
			while (mutaionBit2 == mutaionBit1) {
				mutaionBit2 = rand() % GENE;
			}
			cout << "\tThe second mutation bit: " << setw(2) << mutaionBit2 << endl;

			//Flip the mutation bit (if condition)
			swap(children[c][mutaionBit1], children[c][mutaionBit2]);

		}
		else {
			cout << "\n\tMutation no happen in children " << c + 1<<endl;
		}
	}

	//Print the mutated chromosomes
	cout << "\n";
	cout << "\n\tindex:            0   1   2   3   4   5   6   7   8   9";

	for (int c = 0; c < 2; c++) {
		cout << "\n";
		cout << "\tchildren " << c + 1 << " : " << "\t";
		for (int g = 0; g < GENE; g++) {
			cout << setw(3) << children[c][g] << " ";
		}

	}
	cout << "\n";

	cout << endl;

}

void survivalSelection() {
	//Copy children to the survival chromosome array
	for (int c = 0; c < 2; c++) {
		for (int g = 0; g < GENE; g++) {
			NewChromosome[countNewChromosome][g] = children[c][g];
		}
		//Update array counter
		countNewChromosome++;
	}

	//Print chromosomes in the new population
	cout << "\n\tPrint new chromosome \n";
	for (int c = 0; c < POP_SIZE; c++) {
		cout << "\tchromosome " << c << " :\t";
		for (int g = 0; g < GENE; g++) {
			cout << " " << setw(3) << NewChromosome[c][g];
		}
		cout << "\n";
	}
}

void copyChromosome() {
	cout << "\n\treplace chromosome \n";

	//Copy newChromosome to chromosome
	for (int c = 0; c < POP_SIZE; c++) {
		cout << "\tchromosome " << c << " :\t";
		for (int g = 0; g < GENE; g++) {
			chromosome[c][g] = NewChromosome[c][g];
			cout << " " << setw(3) << chromosome[c][g];
		}
		cout << "\n";
	}
}

void calculateAverageFitness() {
	//1. Declare a variable for totalFitness, initialize to 0
	float totalFitness = 0;

	//2. For every chromosome
	for (int c = 0; c < POP_SIZE; c++) {
		//2.1 Accumulate the fitness into totalFitness
		totalFitness += fitness_value[c];
	}

	//3. Divide the totalFitness with population size
	avgFitness = totalFitness / POP_SIZE;

	//4. Print out the average to the screen
	cout << "\n\tAverage Fitness: " << avgFitness << endl;

	//5. Print out the average to an output file that keep average fitness
	avgFitnessFile << avgFitness << endl;

}

void recordBestFitness() {
	//1. Declare the bestChromosome data structure (global/static)

	//2. For each chromosome

	for (int c = 0; c < POP_SIZE; c++) {

		//2.1. if (fitness current chromosome better than bestFitness){
		if (fitness_value[c] > bestFitness) {
			//2.1.1. bestFitness = fitness for the current chromosome
			bestFitness = fitness_value[c];

			for (int g = 0; g < GENE; g++) {
				//2.1.2. copy the chromosome to bestChromosome
				bestChromosome[g] = chromosome[c][g];
			}

		}

	}

	//3. Print the bestFitness and bestChromosome to the screen
	cout << "\n\tBest Fitness: " << bestFitness << endl;
	cout << "\n\tBest Chromosome:\t";
	for (int g = 0; g < GENE; g++) {
		cout <<" " << setw(3) << bestChromosome[g];
	}
	cout << endl;

	//4. Print the bestFitness and bestChromosome to two separate files
	BestFitnessFile << bestFitness << endl;

	BestChromosomeFile << "\tBest Chromosome:\t";
	for (int g = 0; g < GENE; g++) {
		BestChromosomeFile << " " << setw(3) << bestChromosome[g];
	}
	BestChromosomeFile << endl;
}

void GenerateBestSolutionFile(int Gen) {
	int NumCCTV;
	double SumVisionValue, SumCCTVPrice;
	double f1, f2;
	int allele1;
	string allele2;

	BestSolutionFile << Gen + 1 << ",";
	BestSolutionFile << avgFitness << ",";
	BestSolutionFile << bestFitness << ",,,";

	for (int g = 0; g < GENE; g++) {
		BestSolutionFile << bestChromosome[g] << ",";
	}
	BestSolutionFile << endl;

	//final soluton
	if (Gen+1 == MAXIMUM_GENERATION) {
		BestSolutionFile << "\nFinal Chromosome:,,,,,";
		for (int g = 0; g < GENE; g++) {
			BestSolutionFile << bestChromosome[g] << ",";
		}
		BestSolutionFile << endl;

		//========================================

		NumCCTV = 0;
		f1 = f2 = SumVisionValue = SumCCTVPrice = 0.0;
		initialMap();

		for (int g = 0; g < GENE; g++) {

			if (bestChromosome[g] != "00") {
				NumCCTV++;

				//sample = "3B" , allele1=3, allele2=B 
				if (bestChromosome[g].length() < 3) {
					allele1 = stoi(bestChromosome[g].substr(0, 1));	//Starting position = 0, length = 1 (get num after the position)
					allele2 = bestChromosome[g].substr(1, 1);		//Starting position = 1, length = 1 (get num after the position)
				}
				else {	// need handle when 2 decimal sample = "10A" , allele1=10, allele2=A
					allele1 = stoi(bestChromosome[g].substr(0, 2));	//Starting position = 0, length = 1 (get num after the position)
					allele2 = bestChromosome[g].substr(2, 1);		//Starting position = 1, length = 1 (get num after the position)
				}

				int y = allele1;
				int x = g + 1;	//(1-10)

				//set map vision area
				if (allele2 == "A") {
					//A CCTV
					int r = y - 2;
					int c = x - 2;
					for (int i = r; i < r + 3; i++) {
						for (int j = c; j < c + 3; j++) {
							if ((i >= 0 && i < MAP_SIZE) && (j >= 0 && j < MAP_SIZE)) {
								mapValue[i][j] = 1;
							}
						}
					}
				}
				else if (allele2 == "B") {
					//B CCTV
					int r = y - 3;
					int c = x - 3;
					for (int i = r; i < r + 5; i++) {
						for (int j = c; j < c + 5; j++) {
							if ((i >= 0 && i < MAP_SIZE) && (j >= 0 && j < MAP_SIZE)) {
								mapValue[i][j] = 1;
							}
						}
					}
				}
				else {
					cout << "something went wrong at evaluatebestChromosome";
					exit(0);
				}

				// calculate total CCTV price

				if (allele2 == "A") {
					SumCCTVPrice += CCTV_PRICE[0];
				}
				else if (allele2 == "B") {
					SumCCTVPrice += CCTV_PRICE[1];
				}
				else {
					cout << "something went wrong at calculate total CCTV price";
					exit(0);
				}

			}//end if(bestChromosome[g] != "00")

		}//end one bestChromosome loop

		//display each bestChromosome vision area 2D
		int type = 2;//for Excel
		displayMapValue(type);


		// calculate total vision area
		for (int r = 0; r < MAP_SIZE; r++) {
			for (int c = 0; c < MAP_SIZE; c++) {
				SumVisionValue += mapValue[r][c];
			}
		}

		// calculate fitness 1
		f1 = SumVisionValue / (MAP_SIZE * MAP_SIZE);

		// calculate fitness 2			
		f2 = (max(static_cast<double>(BUDGET), SumCCTVPrice) - SumCCTVPrice) / max(static_cast<double>(BUDGET), SumCCTVPrice);


		// fitness value
		bestFitness = (f1 * 0.5) + (f2 * 0.5);
		BestSolutionFile << "\nNumCCTV:," << NumCCTV << "\nSumVisionValue:," << SumVisionValue << "\nSumCCTVPrice:," << SumCCTVPrice << "\nFitness1:," << f1 << "\nFitness2:," << f2 << "\nTotalFitness:," << bestFitness << endl;


	}


}