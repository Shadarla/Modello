// Modello
using namespace std;
#include <iostream>
#include <vector>
#include <map>
#include <tuple>
#include <math.h>
#include <fstream>
#include <string>
#include <sstream>


//C:\Users\gianc\Documents\GitHub\SSIoT\Sim - n_services_1 - n_devices_25 - n_master_1 - lambda_10.000000 - tot_sim_500 - seed_3 - resource_ctrl_1 - qoe_ctrl_1

double state_probability (int i, int j, int k);
double calcolo_lambda_j();
double trust_model();
void calcolo_denominatore_lambda_k();
double calcolo_lambda_k();
double loss_probability();
double prob_blocco(int index);
double factorial(int n);

int k = 18;
int T = 20;
int nj = 0;
const int Number_of_nodes = 25;
const double lambda = 20;
const int num_classi_di_servizio = 2;
const double mu_j = 1.4; //verifica
const bool type_of_node = 0; //0 per benevolo, 1 per malevolo
const int max_allocab_resources = 2;
const int max_passo = T+7; //esiste una formuletta per calcolarlo
const int num_amici = 13; // supposto uguale per tutti (to check per caso non omogeneo)
const int num_amici_j = 10;
double denominatore_lambda_k = 0; // non varia mai

vector<double> mu;
vector<double> num_amici_per_classe;
vector<double> array_prob_blocco_per_classe;
vector<double> Lambda;


class specifiche_nodo
{
private:
    int id_nodo;
    int num_amici;
    double classe; 
    double mu;
    bool type_of_node;
    double probabilita_feedback_positivo;
     
public:
    vector<double> S; //social factor
    specifiche_nodo();
   // void set_id_nodo(int);
   // int get_id_nodo();
   // void set_num_amici(int);
   // int get_num_amici();
   // void set_probabilità_feedback_positivo(bool);
   // double get_probabilità_feedback_positivo();
  
   

    void set_id_nodo(int id) {
        id_nodo = id;
    }

    int get_id_nodo() {
        return this->id_nodo;
    }

    void set_num_amici(int amici_nodo) {
        num_amici = amici_nodo;
    }

    int get_num_amici() {
        return this->num_amici;
    }


    void set_probabilità_feedback_positivo(bool type_of_node) {
        if (type_of_node == 0) {
            //nodo benevolo
            probabilita_feedback_positivo = 0.9;
        }
        else {
            //nodo malevolo
            probabilita_feedback_positivo = 0.5;
        }
    }

    double get_probabilità_feedback_positivo() {
        return this->probabilita_feedback_positivo;
    }
};

specifiche_nodo::specifiche_nodo() {
    this->id_nodo = 0;
    this->num_amici = 0;
    this->classe = 0;
    this->mu = 0;
    this->type_of_node = 0;
    this->probabilita_feedback_positivo = 0;
}

vector<specifiche_nodo> topologia;


tuple <int, int, int> stato;
map<tuple<int, int, int>, double> mapOfTuple;

int main() {
    int flag_prob_stato = 0;
    int flag_prob_perdita = 0;

    int indice_topologia = 0;
    specifiche_nodo nodo_di_appoggio;
    specifiche_nodo reset;
    double variabile_appoggio = 0;
    string path = "SocialMatrix.txt";
    int contatore_amici = 0;
    string line;
    int social_index = 0;

    ifstream sim_file(path);

    for (indice_topologia = 0; indice_topologia < Number_of_nodes; indice_topologia++) {
        //leggere dato da file e mettere in appoggio
        
        nodo_di_appoggio.set_id_nodo(indice_topologia+1);
             
       //>>

        if (sim_file.is_open()) {

            getline(sim_file, line);
            istringstream iss(line);
                for (social_index = 0; social_index < Number_of_nodes; social_index++)
                {
                    iss >> variabile_appoggio;
                   //cout << variabile_appoggio << '\n';
                    nodo_di_appoggio.S.push_back(variabile_appoggio);
                    if (variabile_appoggio > 0) {
                        contatore_amici++;
                    }
                }       
        }
        nodo_di_appoggio.set_num_amici(contatore_amici);
        contatore_amici = 0;

        // inserire gli altri dati da file
        topologia.push_back(nodo_di_appoggio);
        nodo_di_appoggio = reset;
    }
    sim_file.close();

    if (flag_prob_stato){
    
        stato = make_tuple(k, T, nj);
        double P = 1;
        mapOfTuple[stato] = P;
        calcolo_denominatore_lambda_k();


        cout << "Stato: (" << k << "," << T << "," << nj << ") Probabilita' di stato: " << P << endl;
        // int estract = mapOfTuple[make_tuple(i,j,k)];

        //la prima prob è uno, quindi inizio a calcolare dalla seconda
        nj++;

        for (T = 20; T <= max_passo; T++) {
            for (k = 18; k <= T - 2; k++) {
                for (nj = 0; nj <= max_allocab_resources; nj++) {
                    if (T == 20 && k == 18 && nj == 0)
                        continue;
                    cout << "Stato: (" << k << "," << T << "," << nj << ") ";
                    stato = make_tuple(k, T, nj);
                    P = state_probability(k, T, nj);
                    mapOfTuple[stato] = P;
                    cout << "Probabilita' di stato: " << P << endl;

                }
            }
        }
    }
  
    if (flag_prob_perdita) {

        double prob_perdita = 0;
        prob_perdita = loss_probability();
    }

    return 0;
}

double state_probability(int k, int T, int nj) {
   
    double probabilita = 0;
    double probabilita_feedback_positivo = 0;
    double probabilita_feedback_negativo = 0;
    double lambda_j = calcolo_lambda_j();
    
   
    if (type_of_node == 0) { 
        //nodo benevolo
        probabilita_feedback_positivo = 0.9;
    }
    else { 
        //nodo malevolo
        probabilita_feedback_positivo = 0.5;
    }

    probabilita_feedback_negativo = (1 - probabilita_feedback_positivo);

    // al posto di 0.5 andrà rispettivamente il valore dal map P(kj,Tj,nj-1), P(kj,Tj-1,nj+1), P(kj-1,Tj-1,nj+1)
    //ricordasi di sostuire tutti e tre gli act_res[k] rispettivamente con k-1, e poi tutti k+1 quando si potrà fare

    //controllo sullo stato attraverso map

    double prob_prima_parte = 0;
    if (nj - 1 >= 0){
        prob_prima_parte = mapOfTuple[make_tuple(k, T, nj - 1)] * lambda_j / (lambda_j + (nj - 1) * mu_j);
    }
    else{
        prob_prima_parte = 0;
    }
        
    double prob_seconda_parte = mapOfTuple[make_tuple(k, T-1, nj+1)] * ((nj + 1) * mu_j * probabilita_feedback_negativo) / (lambda_j + (nj + 1) * mu_j);
    double prob_terza_parte = mapOfTuple[make_tuple(k-1, T-1, nj+1)] * ((nj + 1) * mu_j * probabilita_feedback_positivo) / (lambda_j + (nj + 1) * mu_j);

    probabilita = prob_prima_parte + prob_seconda_parte + prob_terza_parte;
    return probabilita;
}

void calcolo_denominatore_lambda_k() {
    int i = 0;
    for (i = 0; i < topologia.size(); i++) {
        denominatore_lambda_k = denominatore_lambda_k + topologia[i].get_probabilità_feedback_positivo() * topologia[i].get_num_amici();
    }
}
 
double calcolo_lambda_k(int indice_amico_di_i) {
    double lambda_k = 1;
    
   //lambda_k = (lambda * probabilita_feedback_positivo * numero amici di k) / denominatore_lambda_k;
    return lambda_k;
}

double calcolo_lambda_j() {
    double lambda_provider = 1;  //lambda_j risultato da ritornare
    vector<double> prob_richiesta_di_i_assegnata_a_j;
    double lambda_ij = (lambda / Number_of_nodes);
    //double prob_soprasoglia = 1;
    //double prob_j_piu_trusted = 1;
    //double prob_j_non_piu_trusted_ma_risorse = 1;
    double trust_check = 0;

    int i;
    for (i = 0; i < num_amici_j-1; i++) {
        trust_check = trust_model(); // dopo teorema delle prob totali sul numero di amici di i piu trusted di j
        prob_richiesta_di_i_assegnata_a_j.push_back(trust_check);
        lambda_provider = lambda_provider + lambda_ij * prob_richiesta_di_i_assegnata_a_j[i];
    }
    

    return lambda_provider;
}

double trust_model() {
    int indice_amico=0;
    double valore_controllo_trust = 0;
    double probabilità_congiunta = 0;
    double prob_amico_piu_trusted_di_j = 1;
    int num_amici_i_piu_trusted = 0; // da verificare
    double lambda_k = 0;
        
    for (indice_amico = 0; indice_amico < num_amici_i_piu_trusted - 1; indice_amico++) {
        if (indice_amico == 0) {
            probabilità_congiunta = 1;
        }
        else {
            lambda_k = calcolo_lambda_k(indice_amico);
            probabilità_congiunta = lambda_k; //prob di blocco del kappesimo amico
        }
        valore_controllo_trust = valore_controllo_trust + (probabilità_congiunta*prob_amico_piu_trusted_di_j);
    }

    return valore_controllo_trust;
}

//CASO NON OMOGENEO
/*double loss_probability() {
    double prob_perdita = 1;
    double traffico_perso = 0;
    int s = 0;

    cout << endl << "Probabilita' di perdita del sistema " << prob_perdita << endl;

    traffico_perso = lambda * prob_perdita;
    cout << endl << "traffico perso dal sistema " << traffico_perso << endl;
    return prob_perdita;
}*/



//CASO OMOGENEO
double loss_probability() {

    double prob_perdita = 1;
    double traffico_perso = 0;
    int s = 0;

    //assegno i valori dei tempi di servizio
    cout << endl << "Stampo valori di mu: " << endl;
    for (s = 0; s < num_classi_di_servizio; s++) {
        if (s == 0) {
            mu.push_back(1.4);
        }
        if (s == 1) {
            mu.push_back(0.7);
        }
        if (s >= 2) {
            mu.push_back(0.025);
        }
        cout << " mu[" << s+1 << "] = " << mu[s] << endl;
    }

    //assegno il numero di amici per ogni classe di servizio

    cout << endl << "Il numero totale di amici per object e': " << num_amici << endl << "Numero di amici per classe: " << endl;
    double proporz_classe1 = 9.3;
    double proporz_classe2 = 3.7;
    double proporz_classe3 = 0;
    double totale_proporzione = proporz_classe1 + proporz_classe2 + proporz_classe3;
    for (s = 0; s < num_classi_di_servizio; s++) {
        if (s == 0) {
            //num_amici_per_classe.push_back(round(num_amici * proporz_classe1 / totale_proporzione));
            num_amici_per_classe.push_back(num_amici * proporz_classe1 / totale_proporzione);
        }
        if (s == 1) {
            //num_amici_per_classe.push_back(round(num_amici * proporz_classe2 / totale_proporzione));
            num_amici_per_classe.push_back(num_amici * proporz_classe2 / totale_proporzione);
        }
        if (s == 2) {
           // num_amici_per_classe.push_back(round(num_amici * proporz_classe3 / totale_proporzione));
            num_amici_per_classe.push_back(num_amici * proporz_classe3 / totale_proporzione);
        }
       //num_amici_per_classe.push_back(num_amici/ num_classi_di_servizio);
        cout << " per classe[" << s+1 << "]: " << num_amici_per_classe[s] << endl;
    }

    Lambda.push_back(lambda / ceil(num_amici_per_classe[0]));
    //Lambda.push_back(lambda / 8);


    cout << endl << "Probabilita' di blocco per servente " << endl;
   // array_prob_blocco_per_classe.push_back(prob_blocco(0));
    //cout << " di classe[1]: " << array_prob_blocco_per_classe[0] << endl;
    for (s = 0; s < num_classi_di_servizio; s++) {
        array_prob_blocco_per_classe.push_back(prob_blocco(s));
        cout << " di classe[" << s+1 << "]: " << array_prob_blocco_per_classe[s] << endl;
        if (s < num_classi_di_servizio - 1)
            Lambda.push_back((lambda * array_prob_blocco_per_classe[s]) / ceil(num_amici_per_classe[s+1])); 
            //Lambda.push_back((lambda * array_prob_blocco_per_classe[s]) / 6);
    }
    
    
    for (s = 0; s < array_prob_blocco_per_classe.size(); s++) {
        prob_perdita = prob_perdita * array_prob_blocco_per_classe[s];
    }

    cout << endl << "Probabilita' di perdita del sistema "<< prob_perdita << endl;
   
    traffico_perso = lambda * prob_perdita;
    cout << endl << "traffico perso dal sistema " << traffico_perso << endl;
    return prob_perdita;
}

double prob_blocco(int index) {
    double P_B_prima_parte = 1;
    double P_B_seconda_parte = 1;
    double P_B_terza_parte = 1;
    double sommatoria = 0;
    double risultato_sommatoria = 0;

    P_B_prima_parte = pow((Lambda[index]/mu[index]),(max_allocab_resources));
    P_B_seconda_parte = 1/factorial(max_allocab_resources);


    int i;
    for (i = 0; i <= (max_allocab_resources); i++) {
        sommatoria = pow((Lambda[index] / mu[index]), i) * (1 / factorial(i));
        risultato_sommatoria = risultato_sommatoria + sommatoria;
    }
    P_B_terza_parte = 1/risultato_sommatoria;

    double P_B = P_B_prima_parte * P_B_seconda_parte * P_B_terza_parte;

    return P_B;
}

double factorial(int n) {
    if (n > 1)
        return n * factorial(n - 1);
    else
        return 1;
}
