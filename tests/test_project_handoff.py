"""Check workflow dependency wiring without importing the container-only SDK."""
import ast
import unittest
from pathlib import Path
from types import SimpleNamespace


class ProjectHandoffTests(unittest.TestCase):
    def test_only_publish_after_all_branches_and_verify_before_cleanup(self):
        tree = ast.parse((Path(__file__).parents[1] / "wf/__init__.py").read_text())
        workflow = next(n for n in tree.body if isinstance(n, ast.FunctionDef))
        workflow.decorator_list = []
        calls = {}

        def task(name):
            def invoke(**kwargs):
                node = SimpleNamespace(name=name, inputs=kwargs)
                calls[name] = node
                if name in ("make_adata", "motifs_task", "publish_archr_project_task"):
                    node.outputs = (object(), object())
                    return node.outputs
                return node
            return invoke

        scope = {"List": list, "Run": object, "Genome": object, "LatchDir": str}
        for node in ast.walk(workflow):
            if isinstance(node, ast.Call) and isinstance(node.func, ast.Name):
                if node.func.id != "LatchDir":
                    scope[node.func.id] = task(node.func.id)
        exec(compile(ast.Module(body=[workflow], type_ignores=[]), "<workflow>", "exec"), scope)
        scope[workflow.name](runs=[], genome="mm10", project_name="example")

        gene = calls["gene_project_task"]
        self.assertIs(calls["motif_coverages_task"].inputs["gene_project_dir"], gene)
        self.assertIs(calls["gene_stats_task"].inputs["gene_project_dir"], gene)
        barrier = calls["complete_results_task"]
        self.assertIs(barrier.inputs["gene_stats_results_dir"], calls["gene_stats_task"])
        self.assertIs(barrier.inputs["motif_results_dir"], calls["motifs_task"].outputs[0])
        publish = calls["publish_archr_project_task"]
        self.assertIs(publish.inputs["results"], barrier)
        self.assertIs(publish.inputs["completed_project"], calls["motifs_task"].outputs[1])
        verify = calls["verify_published_archr_task"]
        self.assertIs(verify.inputs["published_project"], publish.outputs[0])
        self.assertIs(verify.inputs["manifest_hash"], publish.outputs[1])
        self.assertIs(calls["cleanup_checkpoints_task"].inputs["results"], verify)


if __name__ == "__main__":
    unittest.main()
