#include <vector>
#include <iostream>
#include <fstream>
#include <ostream>

using namespace std;

void write_vec_to_file(vector<double> &vect, string filename, ofstream &out, int start, int end);

int main(){
    ofstream out;
    vector<double> vect(100000);
    for (int i = 0; i < 100000; i++){
        vect[i] = i;
        if(i % 10000 == 0 && i > 0){
            cout << "i = "  << i << endl;
            write_vec_to_file(vect, "tstvec.txt", out, i - 10000, i);
            cin.ignore();
        }
    }
    
    write_vec_to_file(vect, "tstvec.txt", out, 0, 100000);
    
    return 0;
}

void write_vec_to_file(vector<double> &vect, string filename, ofstream &out, int start, int end){
    if (start == 0){
        out.open(filename);
    } else {
        out.open(filename, ios::app);
    }
    for (int i = start; i < end; i++){
        out << ' ' << vect[i] << endl;
    }
    out.close();
}
