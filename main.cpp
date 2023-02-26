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


class specifiche_nodo
{
private:
    int id_nodo;
    int num_amici;
    int num_amici_c1;
    int num_amici_c2;
    int classe;
    double mu;
    bool type_of_node;
    double probabilita_feedback_positivo;
    double max_risorse;

public:
    vector<double> S; //social factor
    specifiche_nodo();


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

    void set_num_amici_c1(int amici_nodo_c1) {
        num_amici_c1 = amici_nodo_c1;
    }

    int get_num_amici_c1() {
        return this->num_amici_c1;
    }

    void set_num_amici_c2(int amici_nodo_c2) {
        num_amici_c2 = amici_nodo_c2;
    }

    int get_num_amici_c2() {
        return this->num_amici_c2;
    }


    void set_type(bool tipo) {
        type_of_node = tipo;
    }

    int get_type() {
        return this->type_of_node;
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

    void set_classe(int classe_disp) {
        classe = classe_disp;
    }

    int get_classe() {
        return this->classe;
    }

    void set_mu(int classe_disp) {
        if (classe_disp == 0) {
            //classe piu alta
            mu = 1.4;
        }
        else if (classe_disp == 1) {
            //classe media
            mu = 0.7;
        }
        else if (classe_disp == 2) {
            //classe peggiore
            mu = 0.025;
        }

    }


    double get_mu() {
        return this->mu;
    }

    void set_max_risorse(int classe_disp) {
        if (classe_disp == 0) {
            //classe piu alta
            max_risorse = 3;
        }
        else if (classe_disp == 1) {
            //classe media
            max_risorse = 2;
        }
        else if (classe_disp == 2) {
            //classe peggiore
            max_risorse = 1;
        }

    }


    double get_max_risorse() {
        return this->max_risorse;
    }

};

//costruttore
specifiche_nodo::specifiche_nodo() {
    this->id_nodo = 0;
    this->num_amici = 0;
    this->num_amici_c1 = 0;
    this->num_amici_c2 = 0;
    this->classe = 0;
    this->mu = 0;
    this->type_of_node = 0;
    this->probabilita_feedback_positivo = 0;
    this->max_risorse = 0;
}

//C:\Users\gianc\Documents\GitHub\SSIoT\Sim - n_services_1 - n_devices_25 - n_master_1 - lambda_10.000000 - tot_sim_500 - seed_3 - resource_ctrl_1 - qoe_ctrl_1

double state_probability (int i, int j, int k);
double calcolo_lambda_j();
double trust_model(int id_amico_di_j);
void calcolo_denominatore_lambda_k();
double calcolo_lambda_k(int id_nodo_k);
double loss_probability(double lambda_i);
double prob_blocco(int index);
double factorial(int n);
double prob_binomiale(vector<int> amici_di_i, int id_amico_j);
double prob_binomiale(vector<int> amici_di_i, int id_amico_j, int j);
int proporzionalita_tk(specifiche_nodo km);
double prob_di_blocco_generica(specifiche_nodo km, double lambda_k);

int kj = 18;
int Tj = 20;
int nj = 0;
const int Number_of_nodes = 25;
//const double lambda = 20;
double lambda = 20;
const int num_classi_di_servizio = 2;
const double mu_j = 1.4; //verifica (potrebbe anche non servire piu)
const bool type_of_node = 0; //0 per benevolo, 1 per malevolo (potrebbe anche non servire piu)
const int max_allocab_resources = 2;
const int max_passo = Tj+7; //esiste una formuletta per calcolarlo
const int num_amici = 13; // supposto uguale per tutti (to check per caso non omogeneo e potrebbe non servire piu)
double denominatore_lambda_k = 0; // non varia mai
const int id_nodo_j = 1; //id del nodo da valutare
double lambda_ij = (lambda / Number_of_nodes); //essendo omogenea la distribuzione al momento è uguale per tutti

const double soglia = 0.48;


//vettori per prob blocco
vector<double> mu;
vector<double> num_amici_per_classe;
vector<double> array_prob_blocco_per_classe;
vector<double> Lambda;




vector<specifiche_nodo> topologia;


tuple <int, int, int> stato;
map<tuple<int, int, int>, double> mapOfTuple;

int main() {
    //cosa voglio eseguire?
    int flag_prob_stato = 0;
    int flag_prob_perdita = 1;

    int indice_topologia = 0; //indice per il vettore di classi
    specifiche_nodo nodo_di_appoggio;
    specifiche_nodo reset;
    double variabile_appoggio = 0;
    int variabile_appoggio_per_quel_fallito_di_antonio = 0;
    float appoggio_double=0;
    string path = "SocialMatrix.txt";
    string path2 = "1.UserInfo.txt";
    int contatore_amici = 0;
    string line;
    string line2;
    int social_index = 0;

    //********SETTO LA TOPOLOGIA*************//

    ifstream sim_file(path);
    ifstream sim_file_secondo(path2); //per benevolo malevolo e classe dispositivi
    for (indice_topologia = 0; indice_topologia < Number_of_nodes; indice_topologia++) {
        //leggere dato da file e mettere in appoggio
        
        nodo_di_appoggio.set_id_nodo(indice_topologia+1);
       
        if (sim_file.is_open()) {

            getline(sim_file, line);
            istringstream iss(line);
                for (social_index = 0; social_index < Number_of_nodes; social_index++)
                {
                    iss >> variabile_appoggio;
                   //cout << variabile_appoggio << '\n';
                    nodo_di_appoggio.S.push_back(variabile_appoggio);
                    if (variabile_appoggio > 0) 
                        contatore_amici++;
                }       
        }

        nodo_di_appoggio.set_num_amici(contatore_amici);
        contatore_amici = 0;

        if (sim_file_secondo.is_open()){
            if (indice_topologia == 0) {
                getline(sim_file_secondo, line2);
            }
            getline(sim_file_secondo, line2);
            istringstream iss2(line2);
            for (int ind_ciclo = 0; ind_ciclo < 9; ind_ciclo++) {
                if (ind_ciclo == 6 || ind_ciclo == 7) {
                    iss2 >> appoggio_double;
                    //cout << appoggio_double << endl;
                }
                else {
                    iss2 >> variabile_appoggio_per_quel_fallito_di_antonio;
                    //cout << variabile_appoggio_per_quel_fallito_di_antonio << endl;
                    if (ind_ciclo == 4) {
                        //cout << variabile_appoggio_per_quel_fallito_di_antonio << endl;
                        if (variabile_appoggio_per_quel_fallito_di_antonio == 3)
                            variabile_appoggio_per_quel_fallito_di_antonio = 0;
                        else if (variabile_appoggio_per_quel_fallito_di_antonio == 2)
                            variabile_appoggio_per_quel_fallito_di_antonio = 1;
                        else if (variabile_appoggio_per_quel_fallito_di_antonio == 1)
                            variabile_appoggio_per_quel_fallito_di_antonio = 2;
                        //cout << variabile_appoggio_per_quel_fallito_di_antonio << endl;
                        nodo_di_appoggio.set_classe(variabile_appoggio_per_quel_fallito_di_antonio); 
                    }
                    if (ind_ciclo == 8) {
                        //cout << variabile_appoggio_per_quel_fallito_di_antonio << endl;
                        nodo_di_appoggio.set_type(variabile_appoggio_per_quel_fallito_di_antonio); 
                    }
                }
            }
            nodo_di_appoggio.set_probabilità_feedback_positivo(nodo_di_appoggio.get_type());
            nodo_di_appoggio.set_mu(nodo_di_appoggio.get_classe());
            nodo_di_appoggio.set_max_risorse(nodo_di_appoggio.get_classe());
        }
        topologia.push_back(nodo_di_appoggio);
        nodo_di_appoggio = reset;
    }
    sim_file.close();
    sim_file_secondo.close();

    //************CALCOLO PROBABILITA' DI STATO**************//

    if (flag_prob_stato){
    
        stato = make_tuple(kj, Tj, nj);
        double P = 1;
        mapOfTuple[stato] = P;

        calcolo_denominatore_lambda_k();

        std::cout << "Stato: (" << kj << "," << Tj << "," << nj << ") Probabilita' di stato: " << P << endl;
        // int estract = mapOfTuple[make_tuple(i,j,k)];

        //la prima prob è uno, quindi inizio a calcolare dalla seconda
        nj++;

        for (Tj = 20; Tj <= max_passo; Tj++) {
            for (kj = 18; kj <= Tj - 2; kj++) {
                for (nj = 0; nj <= max_allocab_resources; nj++) {
                    if (Tj == 20 && kj == 18 && nj == 0)
                        continue;
                    std::cout << "Stato: (" << kj << "," << Tj << "," << nj << ") ";
                    stato = make_tuple(kj, Tj, nj);
                    P = state_probability(kj, Tj, nj);
                    mapOfTuple[stato] = P;
                    std::cout << "Probabilita' di stato: " << P << endl;

                }
            }
        }
    }
  
    if (flag_prob_perdita) {

        double prob_perdita = 0;
        double traffico_perso = 0;
   
        int l = 0; //per cicla con lambda diversi

        /*
if (file_risultati.is_open()) {
    //file_risultati << lambda_ij << '\t' << Traffico_perso ;
}
*/
//********SETTO L'OUTPUT DEI RISULTATI*************
     string path_out = "File_output.txt";
     ofstream file_risultati(path_out);
     file_risultati << "LAMBDA" << '\t' << "T_Perso" << '\n';
        for (l = 5; l <= 31; l++) {
            if (file_risultati.is_open()) {
                double lambda_l = ((double)l / (double)Number_of_nodes);
                cout << "Lambda i uguale a " << lambda_l << endl;
            
                //prob_perdita = loss_probability();
                traffico_perso = loss_probability(lambda_l);
                file_risultati << l <<'\t' << traffico_perso << '\n';
            }
        }
     file_risultati.close();
    }

    return 0;
}

double state_probability(int kj, int Tj, int nj) { //forse non serve passare i tre parametri (var globale)
   
    double probabilita = 0; //variabile da ritornare
    double lambda_j = calcolo_lambda_j();
    //double lambda_j = 1;
   
    //controllo sullo stato attraverso map

    double prob_prima_parte = 0;
    if (nj - 1 >= 0){
       // cout << endl << (lambda_j + (nj - 1)) << endl;
        if (lambda_j == 0)
            prob_prima_parte = 0;
        else
            prob_prima_parte = mapOfTuple[make_tuple(kj, Tj, nj - 1)] * lambda_j / (lambda_j + (nj - 1) * topologia[id_nodo_j-1].get_mu());
    }
    else{
        prob_prima_parte = 0;
    }
        
    double prob_seconda_parte = mapOfTuple[make_tuple(kj, Tj-1, nj+1)] * ((nj + 1) * topologia[id_nodo_j-1].get_mu() * (1-topologia[id_nodo_j - 1].get_probabilità_feedback_positivo())) / (lambda_j + (nj + 1) * topologia[id_nodo_j-1].get_mu());
    double prob_terza_parte = mapOfTuple[make_tuple(kj-1, Tj-1, nj+1)] * ((nj + 1) * topologia[id_nodo_j-1].get_mu() * topologia[id_nodo_j - 1].get_probabilità_feedback_positivo()) / (lambda_j + (nj + 1) * topologia[id_nodo_j-1].get_mu());

    //cout << endl << prob_prima_parte << " + " << prob_seconda_parte << " + " << prob_terza_parte << endl;

    probabilita = prob_prima_parte + prob_seconda_parte + prob_terza_parte;
    return probabilita;
}

double calcolo_lambda_j() {
    
    double lambda_provider = 0;  //lambda_j risultato da ritornare
    double prob_richiesta_di_i_assegnata_a_j;

    int i;

    for (i = 0; i < Number_of_nodes; i++) {
        // dopo teorema delle prob totali sul numero di amici di i piu trusted di j
        //if (topologia[id_nodo_j - 1].S[i] > 0) {
        if (topologia[id_nodo_j - 1].S[i] * ((double)kj / (double)Tj) >= soglia) { //se sotto soglia inutile calcolarlo
            // cout << (topologia[id_nodo_j - 1].S[i]) * ((double)kj / (double)Tj) << endl;
            prob_richiesta_di_i_assegnata_a_j = trust_model(i);//passo indice id_amico di j soprasoglia
            //prob_richiesta_di_i_assegnata_a_j = 0.5;//rigo di prova per debuggare
        }
        else
            prob_richiesta_di_i_assegnata_a_j = 0;
        //}
        //else
        //    prob_richiesta_di_i_assegnata_a_j.push_back(0);
        lambda_provider = lambda_provider + (lambda_ij * prob_richiesta_di_i_assegnata_a_j);
    }
    //cout << endl << lambda_provider << endl;
    return lambda_provider;
}


double trust_model(int indice_amico_di_j) {
    double valore_controllo_trust = 0;
    double probabilità_congiunta_k1 = 0;
    //double probabilità_congiunta_k2 = 0;
    double prodotto_prob_congiuta = 0;
    double prob_amico_piu_trusted_di_j = 1;
    double prob_amico_piu_trusted_di_j_parziale = 1;
    //double prob_amico_piu_trusted_di_j2 = 1;
    double prodotto_prob_amico_piu_trusted_di_j = 1;
    int amico_di_i = 0; // da verificare
    double lambda_k = 0;
    //double lambda_k2 = 0;
    int i;
    //int j;
    double T_km;
    int n_index = 0;
    int n_0 = 0;
    double binomiale = 0;
    // int j2;
    specifiche_nodo km;
    vector<int> indici_amici_di_i;

    for (i = 0; i < Number_of_nodes; i++)
    {
        if (topologia[indice_amico_di_j].S[i] > 0 && indice_amico_di_j != (id_nodo_j - 1)) {
            indici_amici_di_i.push_back(i);
        }
    }


    //CONFRONTO AMICI DI I PIU TRUSTED DI J
    //for (i = 0; i < topologia[indice_amico_di_j].get_num_amici()-1; i++) {
    for (i = 0; i < 3; i++) {
        if (i == 0) {
            probabilità_congiunta_k1 = 1;
            //prob_amico_piu_trusted_di_j = 0.5;//DEBUG
            prob_amico_piu_trusted_di_j = prob_binomiale(indici_amici_di_i, indice_amico_di_j); //CONTROLLATA DOVREBBE STARE BENE
            valore_controllo_trust = probabilità_congiunta_k1 * prob_amico_piu_trusted_di_j;
        }
        else if (i == 1) {//RIPRENDERE DA QUI!!
            for (int j = 0; j < indici_amici_di_i.size(); j++) {
                // equivale a dire che il k in questione non ha le risorse disponibili
                //lambda_k = calcolo_lambda_k(indici_amici_di_i[j]);  //calcolo lambda k--QUESTA DA ERRORE
                lambda_k = 1;//DEBUG
                km = topologia[indici_amici_di_i[j]];
                //calcolo T_km
                if (km.get_type() == topologia[id_nodo_j - 1].get_type())
                    T_km = Tj;
                else {
                    T_km = proporzionalita_tk(km);
                }
                probabilità_congiunta_k1 = prob_di_blocco_generica(km, lambda_k); //questo sarà uguale alla prob di blocco di lambda k
                //probabilità_congiunta_k1 = 1;//DEBUG -- CORRETTO
                n_0 = (topologia[id_nodo_j - 1].S[indice_amico_di_j] * ((double)kj / (double)Tj) * T_km) / km.S[indice_amico_di_j];
                //cout << km.S[indice_amico_di_j] << endl;
                //cout << n_0 << endl;
                // DA ERRORE DA CONTROLLARE
                for (n_index = n_0 + 1; n_index <= T_km; n_index++) {
                    prob_amico_piu_trusted_di_j_parziale = prob_binomiale(indici_amici_di_i, indice_amico_di_j, j);//to check termini da passare perche devo passare tutto il vettore tranne indici_amici_di_i[j]
                    //prob_amico_piu_trusted_di_j_parziale = 0.5; //DEBUG
                    binomiale = factorial(T_km) / (factorial(n_index) * factorial(T_km - n_index));//calcolo binomiale tra T_km e n_index;
                    binomiale = binomiale * pow(km.get_probabilità_feedback_positivo(), n_index) * pow(1 - km.get_probabilità_feedback_positivo(), T_km - n_index);
                    prob_amico_piu_trusted_di_j_parziale = prob_amico_piu_trusted_di_j_parziale * binomiale;
                    prob_amico_piu_trusted_di_j = prob_amico_piu_trusted_di_j + prob_amico_piu_trusted_di_j_parziale;
                }
                
                // prob_amico_piu_trusted_di_j = 0.5; //DEBUG
                valore_controllo_trust = valore_controllo_trust + (probabilità_congiunta_k1 * prob_amico_piu_trusted_di_j);
            }
        }
        /*
        caso uguale 2 solo se strettamente necessario
        else if (i == 2) {
            //to check here se posso fermarmi a due
            for (j = 0; j < amici_di_i.size(); j++) {
                lambda_k = calcolo_lambda_k(amici_di_i[j]);
                probabilità_congiunta_k1 = lambda_k; //aggiungere tutta la prob di blocco del kappesimo amico
                for (j2 = 0; j2 < amici_di_i.size(); j2++) {
                    if (amici_di_i[j] != amici_di_i[j2]) {
                        lambda_k2 = calcolo_lambda_k(amici_di_i[j2]);
                        probabilità_congiunta_k2 = lambda_k;
                        prodotto_prob_congiuta = probabilità_congiunta_k1 * probabilità_congiunta_k2;
                        //prob_amico_piu_trusted_di_j
                        //prob_amico_piu_trusted_di_j2

                        //prodotto_prob_amico_piu_trusted_di_j = prob_amico_piu_trusted_di_j * prob_amico_piu_trusted_di_j2;
                        //valore_controllo_trust = valore_controllo_trust + (prodotto_prob_congiuta * prodotto_prob_amico_piu_trusted_di_j);
                    }
                }
            }
        }
        */
    }
    return valore_controllo_trust;
}

void calcolo_denominatore_lambda_k() {
    int i = 0;
    for (i = 0; i < topologia.size(); i++) {
        denominatore_lambda_k = denominatore_lambda_k + topologia[i].get_probabilità_feedback_positivo() * topologia[i].get_num_amici();
    }
}
 
double calcolo_lambda_k(int id_nodo_k) {
    double lambda_k = 1;
    
    lambda_k = (lambda * topologia[id_nodo_k-1].get_probabilità_feedback_positivo() * topologia[id_nodo_k - 1].get_num_amici()) / denominatore_lambda_k;
    return lambda_k;
}

double prob_binomiale(vector<int> indici_amici_di_i,int indice_amico_di_j) {

    double valore_probabilità_parziale = 1;
    double valore_binomiale_sommatoria_totale = 0;
    double binomiale = 0;
    int i = 0;
    int n_index = 0;
    int n0 = 0;
    int id_amico_di_i = 0;
    int T_km = 0;
    specifiche_nodo km;

    for (i = 0; i < indici_amici_di_i.size(); i++) {
        valore_binomiale_sommatoria_totale = 0;
        id_amico_di_i = indici_amici_di_i[i];
        km = topologia[id_amico_di_i];
        if (km.get_type() == topologia[id_nodo_j - 1].get_type())
            T_km = Tj;
        else {
            //to check: formula se non sono dello stesso tipo
            T_km = proporzionalita_tk(km);
        }
        //calcolare n0=(intero piu piccolo Sij*delta j dallo stato*T_km)/Sik;
        n0 = (topologia[id_nodo_j-1].S[indice_amico_di_j] * ((double)kj / (double)Tj) * T_km) / km.S[indice_amico_di_j];

        for (n_index = 0; n_index < n0; n_index++)
        {
            binomiale = factorial(T_km)/ (factorial(n_index) * factorial(T_km- n_index));//calcolo binomiale tra T_km e n_index;
            binomiale = binomiale * pow(km.get_probabilità_feedback_positivo(),n_index) * pow(1 - km.get_probabilità_feedback_positivo(), T_km - n_index);
            valore_binomiale_sommatoria_totale = valore_binomiale_sommatoria_totale + binomiale;
        }
        valore_probabilità_parziale = valore_probabilità_parziale * valore_binomiale_sommatoria_totale;
    }

    return valore_probabilità_parziale;
}

double prob_binomiale(vector<int> indici_amici_di_i, int indice_amico_di_j, int indice_saltare) {

    double valore_probabilità_parziale = 1;
    double valore_binomiale_sommatoria_totale = 0;
    double binomiale = 0;
    int i = 0;
    int n_index = 0;
    int n0 = 0;
    int id_amico_di_i = 0;
    int T_km = 0;
    specifiche_nodo km;

    for (i = 0; i < indici_amici_di_i.size(); i++) {
        if(i != indice_saltare){
            binomiale = 0;
            valore_binomiale_sommatoria_totale = 0;
            id_amico_di_i = indici_amici_di_i[i];
            km = topologia[id_amico_di_i];
            if (km.get_type() == topologia[id_nodo_j - 1].get_type())
                T_km = Tj;
            else {
                //to check: formula se non sono dello stesso tipo
                T_km = proporzionalita_tk(km);
            }
            //calcolare n0=(intero piu piccolo Sij*delta j dallo stato*T_km)/Sik;
            n0 = (topologia[indice_amico_di_j].S[id_nodo_j - 1] * ((double)kj / (double)Tj) * T_km) / km.S[indice_amico_di_j];
            for (n_index = 0; n_index < n0; n_index++)
            {
                binomiale = factorial(T_km) / (factorial(n_index) * factorial(T_km - n_index));//calcolo binomiale tra T_km e n_index;
                binomiale = binomiale * pow(km.get_probabilità_feedback_positivo(), n_index) * pow(1 - km.get_probabilità_feedback_positivo(), T_km - n_index);
                valore_binomiale_sommatoria_totale = valore_binomiale_sommatoria_totale + binomiale;
            }
            valore_probabilità_parziale = valore_probabilità_parziale * valore_binomiale_sommatoria_totale;
        }
    }

    return valore_probabilità_parziale;
}

int proporzionalita_tk(specifiche_nodo km) {
    int T_k = 0;

    //recuperare T_j dallo stato
    T_k = (Tj * km.get_probabilità_feedback_positivo() * km.get_num_amici()) / (topologia[id_nodo_j - 1].get_probabilità_feedback_positivo() * topologia[id_nodo_j - 1].get_num_amici());
    return T_k;

}

double prob_di_blocco_generica(specifiche_nodo km, double lambda_k) {
    double P_B_prima_parte = 1;
    double P_B_seconda_parte = 1;
    double P_B_terza_parte = 1;
    double sommatoria = 0;
    double risultato_sommatoria = 0;

    P_B_prima_parte = pow((lambda_k / km.get_mu()), (km.get_max_risorse()));
    P_B_seconda_parte = 1 / factorial(km.get_max_risorse());


    int i;
    for (i = 1; i <= (km.get_max_risorse()); i++) {  // la i non deve patire da zero?(vedi pdf)
        sommatoria = pow((lambda_k / km.get_mu()), i) * (1 / factorial(i));
        risultato_sommatoria = risultato_sommatoria + sommatoria;
    }
    P_B_terza_parte = 1 / risultato_sommatoria;

    double P_B = P_B_prima_parte * P_B_seconda_parte * P_B_terza_parte;
    //cout << "P_B " <<P_B <<endl; //togli
    return P_B;

}

//CASO NON OMOGENEO
double loss_probability(double lambda_i) {
    double prob_perdita = 1;
    double traffico_perso = 0;
    int s = 0;

    int i = 0;
    int k = 0;
    int j = 0;
    int indice_a = 0;
    int cont_amici_c1 = 0;
    int cont_amici_c2 = 0;
    double Lambda_j_c1 = 0;
    double P_blocco_j_c1 = 0;
    vector<double> Lambda_c1;
    vector<double> P_blocco_c1;

    vector<int> indice_C1;
    vector<int> indice_C2;
    vector<int> indice_C3;


    for (i= 0; i < Number_of_nodes; i++) {

        if (topologia[i].get_classe() == 0) {
            indice_C1.push_back(i);
            //cout << i << endl << endl;
        }
        if (topologia[i].get_classe() == 1) {
            indice_C2.push_back(i);
            //cout << i << endl << endl;
        }
        if (topologia[i].get_classe() == 2) {
            indice_C3.push_back(i);
            //cout << i << endl << endl;
        }
        
        //aggiungo il ciclo per contare gli amici di classe 1 e 2 
        for (indice_a = 0; indice_a < Number_of_nodes; indice_a++) {
            if (topologia[i].S[indice_a] > 0) {
                if (topologia[indice_a].get_classe() == 0) {
                    cont_amici_c1++;
                }
                if (topologia[indice_a].get_classe() == 1) {
                    cont_amici_c2++;
                }
            }
        }

        topologia[i].set_num_amici_c1(cont_amici_c1);
        topologia[i].set_num_amici_c2(cont_amici_c2);
        cont_amici_c1 = 0;
        cont_amici_c2 = 0;
        //cout << "amici c1 del nodo " << i <<": " << topologia[i].get_num_amici_c1() << " \t amici c2 " << topologia[i].get_num_amici_c2() << endl;
    }
    //cout << "DIM indice_c1 " << indice_C1.size() << endl; cout << "DIM indice_c2 " << indice_C2.size() << endl; //prova
    double denom_c1_proporzione = 0;
    int n_1 = 0;
    specifiche_nodo nodo;


   for (j = 0; j < indice_C1.size(); j++) {
       //CALCOLO LAMBDA PER OGNI J
       for (i = 0; i < Number_of_nodes; i++) {
           for (n_1 = 0; n_1 < indice_C1.size(); n_1++) {
               //for (i = 0; i < Number_of_nodes; i++) {
               if (topologia[i].S[indice_C1[n_1]] >= 0) { //per non sommare il contrbuto dise stesso (-1), se è 0 (non amico) somma +0 . somma su tutti gli amici di classe1 grazie al vect ind_C1
                   denom_c1_proporzione = denom_c1_proporzione + (topologia[i].S[indice_C1[n_1]] * topologia[indice_C1[n_1]].get_probabilità_feedback_positivo());
               }
             //}
           }
           if (topologia[indice_C1[j]].S[i] >= 0 && denom_c1_proporzione >0) { //qui faccio filtro solo per i nodi di C1, grzie al vettore Indic_C1, e su se stesso(-1) , 
                Lambda_j_c1 = Lambda_j_c1 + (topologia[indice_C1[j]].S[i] * lambda_i * topologia[indice_C1[j]].get_probabilità_feedback_positivo()/ denom_c1_proporzione );
           }
           denom_c1_proporzione = 0;
       }
       //CALCOLO PROB DI BLOCCO PER OGNI J DATO LAMBDA
       nodo = topologia[indice_C1[j]];
       P_blocco_j_c1 = prob_di_blocco_generica(nodo, Lambda_j_c1);

       //inserimento dei LAMBDA e PROBAB_BLOCCO nei rispettivi vettori "raccoglitori" 
       Lambda_c1.push_back(Lambda_j_c1);
       P_blocco_c1.push_back(P_blocco_j_c1);

       cout <<"Lambda j: " << Lambda_j_c1 << endl << "PB(j): " << P_blocco_j_c1 << endl;
       //cout << denom_c1_proporzione << endl << endl;
       //denom_c1_proporzione = 0;
       Lambda_j_c1 = 0;
       P_blocco_j_c1 = 0;
   }
   
  /* cout << Lambda_c1.size() << endl; //r aggiunte per check
   cout << P_blocco_c1.size() << endl; */
  
   //primo microblocco -- PARTE CLASSE 2
   int indiceJ = 0;
   int indiceK = 0;
   int indice = 0;
   int indice2 = 0;
   int n_2 = 0;
   double Lambda_k_c2 = 0; //riadattare sottto
   vector<double> Lambda_c2;
   double denom_c2_proporzione = 0;
   specifiche_nodo nodo2;
   double P_blocco_k_c2 = 0;
   vector<double> P_blocco_c2;
   double Lambda_k_c2_appoggio = 0;
   double Seconda_parte_appoggio = 0;




   for (indiceK = 0; indiceK < indice_C2.size(); indiceK++) {
       //per il singolo K
       for (indiceJ = 0; indiceJ < Lambda_c1.size(); indiceJ++) {
       
                for (indice = 0; indice < Number_of_nodes; indice++) {
                    for (n_2 = 0; n_2 < indice_C2.size(); n_2++) { // analogo a ciclo con indice n_1
                        if (topologia[indice].S[indice_C2[n_2]] >= 0) { //per non sommare il contrbuto dise stesso (-1), se è 0 (non amico) somma +0 . somma su tutti gli amici di classe2 grazie al vect ind_C2
                            denom_c2_proporzione = denom_c2_proporzione + (topologia[indice].S[indice_C2[n_2]] * topologia[indice_C2[n_2]].get_probabilità_feedback_positivo());
                        }
                    }
                  //(topologia[indice_C2[indiceJ]].S[indice] >= 0)  
                  if (topologia[indice].S[indice_C2[indiceK]] >= 0) {  //qui faccio filtro solo per i nodi di C2, grzie al vettore Indic_C2, e su se stesso(-1) // sotto non deve stare lambdaij, giusto ?
                    Lambda_k_c2_appoggio = Lambda_k_c2_appoggio + (topologia[indice].S[indice_C2[indiceK]] * topologia[indice_C2[indiceK]].get_probabilità_feedback_positivo() / denom_c2_proporzione);
                    //Lambda_k_c2 = Lambda_k_c2 + (Lambda_c1[indiceJ] * P_blocco_c1[indiceJ] * topologia[indice].S[indice_C2[indiceK]] * topologia[indice_C2[indiceK]].get_probabilità_feedback_positivo() / denom_c2_proporzione);
                  } // rimessa la IF
                
                //calcolo la quota parte del traffico delle i  con solo amici di classe2
                  Seconda_parte_appoggio = Seconda_parte_appoggio + lambda_i * topologia[indice].get_num_amici_c2() / max(1, topologia[indice].get_num_amici_c2()) * (1 - min(1, topologia[indice].get_num_amici_c2())) * Lambda_k_c2_appoggio;

                }
                Lambda_k_c2 = Lambda_k_c2 + (Lambda_c1[indiceJ] * P_blocco_c1[indiceJ] * Lambda_k_c2_appoggio) + Seconda_parte_appoggio;
                Lambda_k_c2_appoggio = 0;
                denom_c2_proporzione = 0;
                Seconda_parte_appoggio = 0;
          
       }
        
      //CALCOLO PROB DI BLOCCO PER OGNI K DATO LAMBDA
      nodo2 = topologia[indice_C2[indiceK]];
      P_blocco_k_c2 = prob_di_blocco_generica(nodo2, Lambda_k_c2);
   
      Lambda_c2.push_back(Lambda_k_c2);
      P_blocco_c2.push_back(P_blocco_k_c2);

      cout <<"Stampo per Lambda "<<indiceK <<": " << Lambda_k_c2 << endl;
      cout << "P_blocco per K " << indiceK << ": " << P_blocco_k_c2 << endl;
    
        Lambda_k_c2 = 0;
        P_blocco_k_c2 = 0;
        //cout << Lambda_c2.size() << endl;
        //cout << P_blocco_c2.size() << endl;
   }
   
   //i = 0;
   double T_perso_appoggio = 0;
   double Traffico_perso = 0;

   for (i = 0; i < Lambda_c2.size(); i++) {
       T_perso_appoggio = T_perso_appoggio + (Lambda_c2[i] * P_blocco_c2[i]);
       //cout << endl << "contributo traffico perso " << Lambda_c2[i] * P_blocco_c2[i] << endl;
      // cout << endl << "traffico perso dal sistema " << T_perso_appoggio << endl;
   }
   Traffico_perso = T_perso_appoggio;

   cout << endl << "traffico perso dal sistema " << T_perso_appoggio << endl;

   

   //return prob_perdita; 
   return Traffico_perso;

   /*
    traffico_perso = lambda * prob_perdita;
    std::cout << endl << "traffico perso dal sistema " << traffico_perso << endl;
   */
  


/*
    vector<double> LAMBDA_C1;//ove inserisco i traffici in ingresso a tutti i server di C1
    double Lambda_c1 = 0;
    double Lambda_appoggio = 0;
    vector<double> P_BLOCCO_C1; //ove inserisco le p_blocco di tutti i server di C1
    double P_blocco_c1 = 0;

    

    cout << endl << "Probabilita' di perdita del sistema " << prob_perdita << endl;

    traffico_perso = lambda * prob_perdita;
    cout << endl << "traffico perso dal sistema " << traffico_perso << endl;

    */
    //return prob_perdita;


}

// FINE CASO NON OMOGENEO




 //CASO OMOGENEO
/*
double loss_probability() {

    double prob_perdita = 1;
    double traffico_perso = 0;
    int s = 0;

    //assegno i valori dei tempi di servizio
    std::cout << endl << "Stampo valori di mu: " << endl;
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
        std::cout << " mu[" << s+1 << "] = " << mu[s] << endl;
    }

    //assegno il numero di amici per ogni classe di servizio

    std::cout << endl << "Il numero totale di amici per object e': " << num_amici << endl << "Numero di amici per classe: " << endl;
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
        std::cout << " per classe[" << s+1 << "]: " << num_amici_per_classe[s] << endl;
    }

    Lambda.push_back(lambda / ceil(num_amici_per_classe[0]));
    //Lambda.push_back(lambda / 8);


    std::cout << endl << "Probabilita' di blocco per servente " << endl;
   // array_prob_blocco_per_classe.push_back(prob_blocco(0));
    //cout << " di classe[1]: " << array_prob_blocco_per_classe[0] << endl;
    for (s = 0; s < num_classi_di_servizio; s++) {
        array_prob_blocco_per_classe.push_back(prob_blocco(s));
        std::cout << " di classe[" << s+1 << "]: " << array_prob_blocco_per_classe[s] << endl;
        if (s < num_classi_di_servizio - 1)
            Lambda.push_back((lambda * array_prob_blocco_per_classe[s]) / ceil(num_amici_per_classe[s+1])); 
            //Lambda.push_back((lambda * array_prob_blocco_per_classe[s]) / 6);
    }
    
    
    for (s = 0; s < array_prob_blocco_per_classe.size(); s++) {
        prob_perdita = prob_perdita * array_prob_blocco_per_classe[s];
    }

    std::cout << endl << "Probabilita' di perdita del sistema "<< prob_perdita << endl;
   
    traffico_perso = lambda * prob_perdita;
    std::cout << endl << "traffico perso dal sistema " << traffico_perso << endl;
    return prob_perdita;
}
*/

double prob_blocco(int index) {
    double P_B_prima_parte = 1;
    double P_B_seconda_parte = 1;
    double P_B_terza_parte = 1;
    double sommatoria = 0;
    double risultato_sommatoria = 0;

    P_B_prima_parte = pow((Lambda[index]/mu[index]),(max_allocab_resources));
    P_B_seconda_parte = 1/factorial(max_allocab_resources);


    int i;
    for (i = 1; i <= (max_allocab_resources); i++) {
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
