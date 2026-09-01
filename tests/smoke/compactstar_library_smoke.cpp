// -*- lsst-c++ -*-
/*
 * CompactStar
 * See License file at the top of the source tree.
 *
 * Copyright (c) 2026
 * Mohammadreza Zakeri
 *
 * MIT License — see LICENSE at repo root.
 */

/**
 * @file compactstar_library_smoke.cpp
 * @brief Build/link smoke test for libCompactStar. NOT a scientific test.
 *
 * This program exists to prove that the *build and test mechanism* works:
 *
 *  - CompactStar public headers compile against the installed include paths;
 *  - an executable resolves a genuine out-of-line symbol from libCompactStar.a
 *    (an empty main() linked to a static library can pull no object at all,
 *    which would make the test vacuous);
 *  - the archive and its transitive dependencies link;
 *  - the resulting binary starts and exits cleanly under CTest.
 *
 * It deliberately asserts **no** scientific quantity. There is no stellar mass,
 * moment of inertia, luminosity, EOS value, Hartle result, reaction rate, or
 * reference number anywhere in this file, and none may be added. Scientific
 * validation is roadmap Phase 2A/2B work, governed separately.
 *
 * The chosen symbol is `CompactStar::Core::Prog::GetName()`, defined out-of-line
 * in `CompactStar/Core/src/Prog.cpp`. `Prog` is the non-scientific base class of
 * the library: constructing it performs only member initialisation, it reads no
 * files, opens no network connections, requires no EOS table or data directory,
 * and mutates nothing in the repository.
 */

#include <iostream>
#include <string>

#include "CompactStar/Core/Prog.hpp"

int main()
{
	// Construct through the naming constructor so the name flag is set and
	// GetName() takes its non-logging path.
	const std::string expected = "compactstar_library_smoke";
	CompactStar::Core::Prog prog(expected);

	// Out-of-line call: forces the linker to pull Prog.cpp.o out of the archive.
	const std::string actual = prog.GetName();

	if (actual != expected)
	{
		std::cerr << "compactstar_library_smoke: FAILED — expected \"" << expected
				  << "\", got \"" << actual << "\"\n";
		return 1;
	}

	std::cout << "compactstar_library_smoke: OK (libCompactStar linked and callable)\n";
	return 0;
}
