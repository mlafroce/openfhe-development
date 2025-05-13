//
// Created by matias on 05/11/24.
//


#include "openfhe.h"

using namespace lbcrypto;

constexpr int TEST_MODULUS_SIZE = 5;

NativeVector loadFromBinary(const std::string& path) {
    // Open the file for binary output
    std::ifstream inFile(path, std::ios::binary);
    if (!inFile) {
        throw std::runtime_error("Failed to open file: " + path);
    }

    NativeInteger modulus;
    inFile.read(reinterpret_cast<char*>(&modulus), sizeof(modulus));
    // Write the size of the vector to the file
    size_t size;
    inFile.read(reinterpret_cast<char*>(&size), sizeof(size));

    NativeVector vec(size, modulus);
    // Write the contents of the vector to the file
    inFile.read(reinterpret_cast<char*>(&vec[0]), size * sizeof(NativeInteger));

    // Check for write errors
    if (!inFile) {
        throw std::runtime_error("Failed to write to file: " + path);
    }
    return vec;
}

void encryptVector() {
    // Sample Program: Step 1: Set CryptoContext
    CCParams<CryptoContextBFVRNS> parameters;
    parameters.SetPlaintextModulus(65537);
    parameters.SetMultiplicativeDepth(2);

    CryptoContext<DCRTPoly> cryptoContext = GenCryptoContext(parameters);
    // Enable features that you wish to use
    cryptoContext->Enable(PKE);
    cryptoContext->Enable(KEYSWITCH);
    cryptoContext->Enable(LEVELEDSHE);

    // Sample Program: Step 2: Key Generation

    // Initialize Public Key Containers
    KeyPair<DCRTPoly> keyPair;

    // Generate a public/private key pair
    keyPair = cryptoContext->KeyGen();

    // Generate the relinearization key
    cryptoContext->EvalMultKeyGen(keyPair.secretKey);

    // Generate the rotation evaluation keys
    cryptoContext->EvalRotateKeyGen(keyPair.secretKey, {1, 2, -1, -2});

    // Sample Program: Step 3: Encryption

    // First plaintext vector is encoded
    std::vector<int64_t> vectorOfInts1(0x100);
    for (size_t i = 0; i<vectorOfInts1.size(); ++i)
    {
        vectorOfInts1[i] = i * 2;
    }
    Plaintext plaintext1               = cryptoContext->MakePackedPlaintext(vectorOfInts1);
    std::cout << "Encrypted: " << plaintext1 << std::endl;
    Plaintext plaintextDec;
    auto ciphertext1 = cryptoContext->Encrypt(keyPair.publicKey, plaintext1);
    cryptoContext->Decrypt(keyPair.secretKey, ciphertext1, &plaintextDec);
    std::cout << "Decrypted: " << plaintextDec << std::endl;
}

void runCRT() {
    size_t m = 16;
    size_t n = 8;
    NativeInteger modulusQ(LastPrime<NativeInteger>(TEST_MODULUS_SIZE, m));
    std::cout << "Modulus: " << modulusQ << std::endl;
    // con m = 4 q = 29 (último primo de 5 bits)
    // ya que (q-1) / m es entero (28/4 = 8)
    // con m = 16 q = 17
    NativeInteger rootOfUnity = RootOfUnity(m, modulusQ);
    std::cout << "rootOfUnity: " << rootOfUnity << std::endl;

    //DiscreteUniformGeneratorImpl<NativeVector> dug;
    NativeVector vec(n, modulusQ);
    for (size_t i = 0; i<n; ++i)
    {
        vec[i] = 1;
    }

    std::cout << "Generated vector: " << vec << std::endl;

    ChineseRemainderTransformFTT<NativeVector> crtFTT;
    crtFTT.PreCompute(rootOfUnity, m, modulusQ);

    //crtFTT.InverseTransformFromBitReverseInPlace(rootOfUnity, m, &vec);
    crtFTT.ForwardTransformToBitReverseInPlace(rootOfUnity, m, &vec);

    intnat::saveToBinary(vec, std::string("output.bin"));
    NativeVector vec2 = loadFromBinary(std::string("output.bin"));
    
    std::cout << "Transformed vector: " << vec << std::endl;
    std::cout << "Loaded vector: " << vec2 << std::endl;
}

int main() {
    //runCRT();
    encryptVector();
    return 0;
}