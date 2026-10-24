import sys
import unittest
from pathlib import Path
from unittest.mock import patch
import io

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from experiments import run_all_experiments
from experiments import case_b_symbolic
from experiments import verify_lift_proposition
from experiments import verify_theorem2_random


class TestExperiments(unittest.TestCase):
    @patch('sys.stdout', new_callable=io.StringIO)
    def test_run_all_experiments(self, mock_stdout):
        run_all_experiments.main()
        self.assertIn("Table 1", mock_stdout.getvalue())

    @patch('sys.stdout', new_callable=io.StringIO)
    def test_case_b_symbolic(self, mock_stdout):
        case_b_symbolic.main()
        self.assertIn("PASSED", mock_stdout.getvalue())

    @patch('sys.stdout', new_callable=io.StringIO)
    def test_verify_theorem2_random(self, mock_stdout):
        import experiments.verify_theorem2_random
        # Override samples to make it fast for testing
        experiments.verify_theorem2_random.SAMPLES = 50
        result = experiments.verify_theorem2_random.main()
        self.assertEqual(result, 0)
        self.assertIn("RESULT: PASS", mock_stdout.getvalue())

if __name__ == '__main__':
    unittest.main()
