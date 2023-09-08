#!/usr/bin/env python3

import unittest

from sd2ba.sd2ba_functions.align import start_pos_in_trimmed_sequence


class TestAlign(unittest.TestCase):

    def test_align(self):

        complete = "MFKKLLHSLFAGLTFMAAVAAVPAHAQEADAQTTVKTAVDDVLATIKGDSDLRSGNMQKVFQLVDQKIVPRADFKRTTQIAMGRFWSQATPEQQQQIQDGFKTLLVRTYAGALANVRNQTVSYKPFRAAADDTDVVVRSTVNNNGEPVALDYRMEKSANGWKVYDINISGLWLSETYKNQFAEVISKRGGVSGLVQFLDERNAQLAKGPAK"

        trimmed = "-TAVDDVLATIKGDSDLRSGNMQKVFQLVDQKIVPRADFKRTTQIAMGRFWSQATPEQQQQIQDGFKTLLVRTYAGALANVR-NQTVSYKPFRAAADD-TDVVVRSTVNNNG--EPVALDYRMEKSANGWKVYDINISGLWLSETYKNQFAEVISKRgGVSGLVQFLD-"

        # Trimmed sequence starts at position 35 in the complete sequence
        self.assertEqual(start_pos_in_trimmed_sequence(complete, trimmed), 35)

if __name__ == "__main__":
    unittest.main()
