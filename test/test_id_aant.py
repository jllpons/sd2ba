#!/usr/bin/env python3

import unittest

from sd2ba.id_aant import FetchData

class TestFetchData(unittest.TestCase):

    def setUp(self):
        self.fetchdata = FetchData()

    def test_get_uniprot_entry_data(self):

        data = self.fetchdata.get_uniprot_entry_data('Q8XV73')

        sequence = "MFKKLLHSLFAGLTFMAAVAAVPAHAQEADAQTTVKTAVDDVLATIKGDSDLRSGNMQKVFQLVDQKIVPRADFKRTTQIAMGRFWSQATPEQQQQIQDGFKTLLVRTYAGALANVRNQTVSYKPFRAAADDTDVVVRSTVNNNGEPVALDYRMEKSANGWKVYDINISGLWLSETYKNQFAEVISKRGGVSGLVQFLDERNAQLAKGPAK"

        self.assertEqual(data.data["uniprot_accession"], 'Q8XV73')
        self.assertEqual(data.data["ena_accession"], 'CAD16665')
        self.assertEqual(data.data["uniprot_aa_sequence"], sequence)

