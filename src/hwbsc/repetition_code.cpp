#include <iostream>
#include <cmath>
#include <vector>
#include <numeric>
#include <cstdlib>
#include <ctime>
using namespace std;

vector<vector<int>> repetition_code_encoder(int repetition_code_length, vector<int> binary_string) {
    vector<vector<int>> encoding_output(binary_string.size(), vector<int>(repetition_code_length,0));
    
    for (int i = 0; i < binary_string.size(); i++) {
        for (int j = 0; j < repetition_code_length; j++) {
            encoding_output[i][j] = binary_string[i];
        }

    } 
    return encoding_output;
}

vector<vector<int>> apply_error(vector<vector<int>> encoded_string, float probability) {
    //function that takes the encoded string and applies a random binary error to it 

    vector<vector<int>> error(encoded_string.size(), vector<int>(encoded_string[0].size(),0));
    vector<vector<int>> flawed_output = encoded_string;
    // Seed the random number generator with the current time
    srand(time(0));
    
    for (int i = 0; i < encoded_string.size(); i++) {
        for (int j = 0; j < encoded_string[0].size(); j++) {
            // Generate a random number between 0 and 1
            double random_number = (double) rand() / RAND_MAX;
            if (random_number < probability) {
                error[i][j] = 1;
            }
            flawed_output[i][j] ^= error[i][j];
        }
    }   
    return flawed_output;
}


vector<vector<int>> repetition_code_majority_vote_decoder(vector<vector<int>> binary_encoding)  {
    vector<vector<int>> decoded_output = binary_encoding;
    for (int i= 0; i < binary_encoding.size(); i++) {
        if (int sum = accumulate(binary_encoding[i].begin(), binary_encoding[i].end(), 0) > binary_encoding[0].size()/2) {
            for (int j = 0; j < binary_encoding[0].size() ; j++) {
                decoded_output[i][j] = 1;
            }
        }
        else {
            for (int j = 0; j < binary_encoding[0].size() ; j++) {
                decoded_output[i][j] = 0;
            }
        }
    }
    return decoded_output;
}

int logical_error_counter(vector<int> binary_string, vector<vector<int>> decoded_output) {
    int counter = 0;
    for (int i = 0; i < binary_string.size(); i++)  {
        if (binary_string[i] != decoded_output[i][0]) {
            counter += 1;
        }
    }
    return counter;
}


int main() {
    int repetition_code_length = 3;
    vector<int> bitstring = {1, 0, 1, 0, 1, 0, 1};
    
    vector<vector<int>> encoded = repetition_code_encoder(repetition_code_length, bitstring);
    cout << "normal" << endl;
    for (int i = 0; i < encoded.size(); i++) {
        for (int j = 0; j < encoded[i].size(); j++) {
            cout <<  encoded[i][j] << " ";
        }
        cout << endl;
    }

    vector<vector<int>> flawed = apply_error(encoded , 0.10);
    cout << "flawed" << endl;
    for (int i = 0; i < flawed.size(); i++) {
        for (int j = 0; j < flawed[i].size(); j++) {
            cout <<  flawed[i][j] << " ";
        }
        cout << endl;
    }

    vector<vector<int>> prediction = repetition_code_majority_vote_decoder(flawed);
    cout << "fixed" << endl;
    for (int i = 0; i < prediction.size(); i++) {
        for (int j = 0; j < prediction[i].size(); j++) {
            cout <<  prediction[i][j] << " ";
        }
        cout << endl;
    }

    int counter = logical_error_counter(bitstring, prediction);
    cout << counter << endl;

    return 0;

}