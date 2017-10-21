from motif_remover import *

class TestesDeMontagemDeSequencia(unittest.TestCase):
    def test_calcular_frame_0(self):
        frame = frame_calc("aaa")
        self.assertEqual(frame, 0)

    def test_calcular_frame_1(self):
        frame = frame_calc("aaaa")
        self.assertEqual(frame, 1)

    def test_calcular_frame_2(self):
        frame = frame_calc("aaaaa")
        self.assertEqual(frame, 2)


class TesteDeCalculoDeFrame(unittest.TestCase):
    def setUp(self):
        self.seq0 = "aaaxxxxxxaaaaxxxxxxaaaaxxxxxx"
        self.seq1 = "aaaxxxxxxaaaaaxxxxxx"
        self.seq2 = "actgxxxxxxactgaxxxxxxactg"
        self.seq3 = "xxxxxxactgaxxxxxxactg"

    def test_seq0_substitutions(self):
        output = motif_frame_counter(self.seq0)
        expected_output = [0, 1, 2]
        self.assertEqual(output, expected_output)

    def test_seq1_substitutions(self):
        output = motif_frame_counter(self.seq1)
        expected_output = [0, 2]
        self.assertEqual(output, expected_output)

    def test_seq2_substitutions(self):
        output = motif_frame_counter(self.seq2)
        expected_output = [1, 0]
        self.assertEqual(output, expected_output)

    def test_seq3_substitutions(self):
        output = motif_frame_counter(self.seq3)
        expected_output = [0, 2]
        self.assertEqual(output, expected_output)


class TesteDeReplaceDeMotivo(unittest.TestCase):
    def setUp(self):
        self.seq0 = "aaaxxxxxxaaaaxxxxxxaaaaxxxxxx"
        self.seq1 = "aaaxxxxxxaaaaaxxxxxx"
        self.seq2 = "actgxxxxxxactgaxxxxxxactg"
        self.seq3 = "xxxxxxactgaxxxxxxactg"

    def test_seq0_replace(self):
        output = motif_replacer(self.seq0, "0", "1", "2")
        expected_output = "aaa0aaaa1aaaa2"
        self.assertEqual(output, expected_output)

    def test_seq1_replace(self):
        output = motif_replacer(self.seq1, "0", "1", "2")
        expected_output = "aaa0aaaaa2"
        self.assertEqual(output, expected_output)

    def test_seq2_replace(self):
        output = motif_replacer(self.seq2, "0", "1", "2")
        expected_output = "actg1actga0actg"
        self.assertEqual(output, expected_output)

    def test_seq3_replace(self):
        output = motif_replacer(self.seq3, "0", "1", "2")
        expected_output = "0actga2actg"
        self.assertEqual(output, expected_output)

if __name__ == '__main__':
    unittest.main()