#include "openfhe.h"

extern "C" NativeInteger mulModWrapper(NativeInteger a, NativeInteger b, NativeInteger mod) {
	return a.ModMul(b, mod);
}

int main(int argc, char** argv) {
	NativeInteger a(0x1000200030004000);
	NativeInteger b(0x2000000000000001);
	NativeInteger mod(0x10000);
	NativeInteger c{mulModWrapper(a, b, mod)};
	std::cout << std::hex << "c " << c << std::endl;
}
