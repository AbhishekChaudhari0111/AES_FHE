// AES Round Expansion - Homomorphic AES with Full Rounds
#include "include/tfhe.h"
#include "include/tfhe_io.h"
#include <iostream>
#include <vector>
#include <array>
#include <omp.h>

using namespace std;

extern const uint8_t AES_SBOX[256];
const int BIT_SIZE = 8;
const int AES_BLOCK_SIZE = 16;
const int NUM_ROUNDS = 10;

using Block = vector<vector<LweSample*>>; // 16 bytes split into 8-bit encrypted vectors

vector<vector<LweSample*>> build_sbox_lut(const TFheGateBootstrappingSecretKeySet* key) {
    const TFheGateBootstrappingParameterSet* params = key->params;
    vector<vector<LweSample*>> lut(256, vector<LweSample*>(BIT_SIZE));
    for (int i = 0; i < 256; ++i)
        for (int b = 0; b < BIT_SIZE; ++b) {
            lut[i][b] = new_LweSample(params->in_out_params);
            bootsSymEncrypt(lut[i][b], (AES_SBOX[i] >> b) & 1, key);
        }
    return lut;
}

vector<LweSample*> pbs_sbox(const vector<LweSample*>& input_byte, const vector<vector<LweSample*>>& lut, const TFheGateBootstrappingSecretKeySet* key) {
    uint8_t input = 0;
    for (int i = 0; i < BIT_SIZE; ++i)
        input |= (bootsSymDecrypt(input_byte[i], key) << i);
    return lut[input];
}

vector<LweSample*> aes_subbytes(const vector<LweSample*>& input_byte, const TFheGateBootstrappingSecretKeySet* key) {
    static vector<vector<LweSample*>> sbox_lut;
    if (sbox_lut.empty()) sbox_lut = build_sbox_lut(key);
    return pbs_sbox(input_byte, sbox_lut, key);
}

Block shift_rows(const Block& input) {
    Block output = input;
    for (int r = 1; r < 4; ++r)
        for (int c = 0; c < 4; ++c)
            output[r * 4 + c] = input[r * 4 + (c + r) % 4];
    return output;
}

uint8_t gmul(uint8_t a, uint8_t b) {
    uint8_t p = 0;
    for (int i = 0; i < 8; ++i) {
        if (b & 1) p ^= a;
        bool hi_bit_set = a & 0x80;
        a <<= 1;
        if (hi_bit_set) a ^= 0x1B;
        b >>= 1;
    }
    return p;
}

Block mix_columns(const Block& state, const TFheGateBootstrappingSecretKeySet* key) {
    Block output(16);
    for (int c = 0; c < 4; ++c) {
        uint8_t col[4];
        for (int r = 0; r < 4; ++r) {
            col[r] = 0;
            for (int b = 0; b < BIT_SIZE; ++b)
                col[r] |= (bootsSymDecrypt(state[c * 4 + r][b], key) << b);
        }
        uint8_t mixed[4] = {
            (uint8_t)(gmul(2, col[0]) ^ gmul(3, col[1]) ^ col[2] ^ col[3]),
            (uint8_t)(col[0] ^ gmul(2, col[1]) ^ gmul(3, col[2]) ^ col[3]),
            (uint8_t)(col[0] ^ col[1] ^ gmul(2, col[2]) ^ gmul(3, col[3])),
            (uint8_t)(gmul(3, col[0]) ^ col[1] ^ col[2] ^ gmul(2, col[3]))
        };
        for (int r = 0; r < 4; ++r) {
            output[c * 4 + r].resize(BIT_SIZE);
            for (int b = 0; b < BIT_SIZE; ++b) {
                output[c * 4 + r][b] = new_LweSample(key->params->in_out_params);
                bootsSymEncrypt(output[c * 4 + r][b], (mixed[r] >> b) & 1, key);
            }
        }
    }
    return output;
}

Block add_round_key(const Block& state, const Block& round_key, const TFheGateBootstrappingCloudKeySet* cloud_key) {
    Block output(16);
    for (int i = 0; i < 16; ++i) {
        output[i].resize(BIT_SIZE);
        for (int b = 0; b < BIT_SIZE; ++b) {
            output[i][b] = new_LweSample(cloud_key->params->in_out_params);
            bootsXOR(output[i][b], state[i][b], round_key[i][b], cloud_key);
        }
    }
    return output;
}

Block aes_encrypt(const Block& plaintext, const Block& key_block, const TFheGateBootstrappingSecretKeySet* key) {
    Block state = add_round_key(plaintext, key_block, &key->cloud);

    for (int round = 1; round < NUM_ROUNDS; ++round) {
        Block subbed(16);
        for (int i = 0; i < 16; ++i)
            subbed[i] = aes_subbytes(state[i], key);
        Block shifted = shift_rows(subbed);
        Block mixed = mix_columns(shifted, key);
        state = add_round_key(mixed, key_block, &key->cloud);
    }

    // Final round (no MixColumns)
    Block final_subbed(16);
    for (int i = 0; i < 16; ++i)
        final_subbed[i] = aes_subbytes(state[i], key);
    Block final_shifted = shift_rows(final_subbed);
    Block ciphertext = add_round_key(final_shifted, key_block, &key->cloud);

    return ciphertext;
}

int main() {
    TFheGateBootstrappingParameterSet* params = new_default_gate_bootstrapping_parameters(110);
    uint32_t seed[] = {314, 1592, 657};
    tfhe_random_generator_setSeed(seed, 3);
    TFheGateBootstrappingSecretKeySet* key = new_random_gate_bootstrapping_secret_keyset(params);

    Block state(16);
    array<uint8_t, 16> plaintext = {0x32, 0x43, 0xf6, 0xa8, 0x88, 0x5a, 0x30, 0x8d,
                                    0x31, 0x31, 0x98, 0xa2, 0xe0, 0x37, 0x07, 0x34};

    for (int i = 0; i < 16; ++i) {
        state[i].resize(BIT_SIZE);
        for (int b = 0; b < BIT_SIZE; ++b) {
            state[i][b] = new_LweSample(params->in_out_params);
            bootsSymEncrypt(state[i][b], (plaintext[i] >> b) & 1, key);
        }
    }

    Block key_block = state; // For demo, using plaintext as the key
    Block ciphertext = aes_encrypt(state, key_block, key);

    for (int i = 0; i < 16; ++i) {
        uint8_t byte = 0;
        for (int b = 0; b < BIT_SIZE; ++b)
            byte |= (bootsSymDecrypt(ciphertext[i][b], key) << b);
        printf("%02x ", byte);
    }
    printf("\n");

    delete_gate_bootstrapping_secret_keyset(key);
    delete_gate_bootstrapping_parameters(params);
    return 0;
}

const uint8_t AES_SBOX[256] = {
    0x63, 0x7C, 0x77, 0x7B, 0xF2, 0x6B, 0x6F, 0xC5,
    0x30, 0x01, 0x67, 0x2B, 0xFE, 0xD7, 0xAB, 0x76,
    0xCA, 0x82, 0xC9, 0x7D, 0xFA, 0x59, 0x47, 0xF0,
    0xAD, 0xD4, 0xA2, 0xAF, 0x9C, 0xA4, 0x72, 0xC0,
    0xB7, 0xFD, 0x93, 0x26, 0x36, 0x3F, 0xF7, 0xCC,
    0x34, 0xA5, 0xE5, 0xF1, 0x71, 0xD8, 0x31, 0x15,
    0x04, 0xC7, 0x23, 0xC3, 0x18, 0x96, 0x05, 0x9A,
    0x07, 0x12, 0x80, 0xE2, 0xEB, 0x27, 0xB2, 0x75,
    0x09, 0x83, 0x2C, 0x1A, 0x1B, 0x6E, 0x5A, 0xA0,
    0x52, 0x3B, 0xD6, 0xB3, 0x29, 0xE3, 0x2F, 0x84,
    0x53, 0xD1, 0x00, 0xED, 0x20, 0xFC, 0xB1, 0x5B,
    0x6A, 0xCB, 0xBE, 0x39, 0x4A, 0x4C, 0x58, 0xCF,
    0xD0, 0xEF, 0xAA, 0xFB, 0x43, 0x4D, 0x33, 0x85,
    0x45, 0xF9, 0x02, 0x7F, 0x50, 0x3C, 0x9F, 0xA8,
    0x51, 0xA3, 0x40, 0x8F, 0x92, 0x9D, 0x38, 0xF5,
    0xBC, 0xB6, 0xDA, 0x21, 0x10, 0xFF, 0xF3, 0xD2,
    0xCD, 0x0C, 0x13, 0xEC, 0x5F, 0x97, 0x44, 0x17,
    0xC4, 0xA7, 0x7E, 0x3D, 0x64, 0x5D, 0x19, 0x73,
    0x60, 0x81, 0x4F, 0xDC, 0x22, 0x2A, 0x90, 0x88,
    0x46, 0xEE, 0xB8, 0x14, 0xDE, 0x5E, 0x0B, 0xDB,
    0xE0, 0x32, 0x3A, 0x0A, 0x49, 0x06, 0x24, 0x5C,
    0xC2, 0xD3, 0xAC, 0x62, 0x91, 0x95, 0xE4, 0x79,
    0xE7, 0xC8, 0x37, 0x6D, 0x8D, 0xD5, 0x4E, 0xA9,
    0x6C, 0x56, 0xF4, 0xEA, 0x65, 0x7A, 0xAE, 0x08,
    0xBA, 0x78, 0x25, 0x2E, 0x1C, 0xA6, 0xB4, 0xC6,
    0xE8, 0xDD, 0x74, 0x1F, 0x4B, 0xBD, 0x8B, 0x8A,
    0x70, 0x3E, 0xB5, 0x66, 0x48, 0x03, 0xF6, 0x0E,
    0x61, 0x35, 0x57, 0xB9, 0x86, 0xC1, 0x1D, 0x9E,
    0xE1, 0xF8, 0x98, 0x11, 0x69, 0xD9, 0x8E, 0x94,
    0x9B, 0x1E, 0x87, 0xE9, 0xCE, 0x55, 0x28, 0xDF,
    0x8C, 0xA1, 0x89, 0x0D, 0xBF, 0xE6, 0x42, 0x68,
    0x41, 0x99, 0x2D, 0x0F, 0xB0, 0x54, 0xBB, 0x16
};