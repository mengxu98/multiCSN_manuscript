import csv
import gzip
import shutil
import sys
import tempfile
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "functions"))
from runtime import REPO
from file_io import publish_directory

SPECIFICATIONS = [
    ('results/gse194122/cell_budget/module_coupling.tsv', 'erythroid_cell_budget/module_coupling.tsv', 'Figure 5e'),
    ('results/gse194122/cell_budget/module_state_profiles.tsv.gz', 'erythroid_cell_budget/module_state_profiles.tsv.gz', 'Figure 5e'),
    ('results/gse194122/cell_budget/focal_GO_summary.tsv', 'erythroid_cell_budget/focal_GO_summary.tsv', 'Methods: Erythroid network stability and program analysis'),
    ('results/gse194122/cell_budget/focal_GO_all_tests.tsv.gz', 'erythroid_cell_budget/focal_GO_all_tests.tsv.gz', 'Methods: Erythroid network stability and program analysis'),
    ('results/gse194122/cases/HMGA1_gene_recurrence.tsv', 'erythroid_application_validation/HMGA1_gene_recurrence.tsv', 'Figure 5i; Figure S15a; Methods: Erythroid network stability and program analysis'),
    ('results/gse194122/cases/HMGA1_GO_all_tests.tsv.gz', 'erythroid_application_validation/HMGA1_GO_all_tests.tsv.gz', 'Figure 5i; Methods: Erythroid network stability and program analysis'),
    ('results/gse194122/cases/HMGA1_size_only_reference.tsv', 'erythroid_application_validation/HMGA1_size_only_reference.tsv', 'Methods: Erythroid network stability and program analysis'),
    ('results/gse194122/cases/GATA1_all_reference_chain_tests.tsv.gz', 'erythroid_application_validation/GATA1_all_reference_chain_tests.tsv.gz', 'Table S6; Methods'),
    ('results/gse194122/cases/GATA1_all_timepoint_tests.tsv.gz', 'erythroid_application_validation/GATA1_all_timepoint_tests.tsv.gz', 'Table S6; Methods'),
    ('results/gse194122/cases/focal_TF_raw_expression.tsv', 'erythroid_functional_programs/focal_TF_raw_expression.tsv', 'Figure 5i,j; Figure S15c'),
    ('results/gse194122/programs/direct_module_GO_all_tests.tsv.gz', 'erythroid_regulatory_programs/direct_module_GO_all_tests.tsv.gz', 'Figure 5j; Figure S15c; Methods: Erythroid network stability and program analysis'),
    ('results/gse194122/programs/chain_module_GO_all_tests.tsv.gz', 'erythroid_regulatory_programs/chain_module_GO_all_tests.tsv.gz', 'Figure S15d; Methods: Erythroid network stability and program analysis'),
    ('results/gse194122/programs/all_module_RNA_summary.tsv', 'erythroid_regulatory_programs/all_module_RNA_summary.tsv', 'Methods: Erythroid network stability and program analysis'),
    ('results/gse194122/programs/all_module_RNA_contrasts.tsv', 'erythroid_regulatory_programs/all_module_RNA_contrasts.tsv', 'Methods: Erythroid network stability and program analysis'),
    ('results/gse194122/programs/HMGA1_BNC2_background_GO_all_tests.tsv.gz', 'erythroid_regulatory_programs/HMGA1_BNC2_background_GO_all_tests.tsv.gz', 'Methods: Erythroid network stability and program analysis'),
    ('results/gse194122/cases/all_TF_state_support.tsv.gz', 'erythroid_link_recurrence/all_TF_state_support.tsv.gz', 'Figure 5f; Figure S14; Methods: Erythroid network stability and program analysis'),
    ('results/gse194122/cases/focal_TF_transition_sizes.tsv', 'erythroid_link_recurrence/focal_TF_transition_sizes.tsv', 'Figure 5i; Figure S15d; Methods: Erythroid network stability and program analysis'),
    ('results/gse194122/cases/all_TF_expression_by_state.tsv', 'erythroid_link_recurrence/all_TF_expression_by_state.tsv', 'Methods: Erythroid network stability and program analysis'),
    ('results/gse194122/cases/all_TF_paired_expression_changes.tsv', 'erythroid_link_recurrence/all_TF_paired_expression_changes.tsv', 'Methods: Erythroid network stability and program analysis'),
    ('results/gse162170/locus_display_ordering.tsv', 'locus_display_ordering.tsv', 'Table S4'),
]

PRIMARY_TERM_SELECTIONS = {
    'results/gse194122/cell_budget/focal_GO_all_tests.tsv.gz': (
        'results/gse194122/cell_budget/displayed_GO_tests.tsv',
        ['seed', 'regulator', 'transition_order', 'topology_class', 'go_id']),
    'results/gse194122/cases/HMGA1_GO_all_tests.tsv.gz': (
        'results/gse194122/cases/HMGA1_displayed_GO_tests.tsv',
        ['background', 'set', 'go_id']),
    'results/gse194122/programs/HMGA1_BNC2_background_GO_all_tests.tsv.gz': (
        'results/gse194122/programs/HMGA1_BNC2_original_terms.tsv',
        ['regulator', 'transition_order', 'event', 'background', 'go_id']),
}


def export_go_tests(source, selected, destination, keys):
    with selected.open() as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        selected_fields = reader.fieldnames
        selected_rows = {tuple(row[key] for key in keys): row for row in reader}
    with gzip.open(source, "rt") as handle, gzip.open(destination, "wt") as output:
        reader = csv.DictReader(handle, delimiter="\t")
        fields = reader.fieldnames + [field for field in selected_fields if field not in reader.fieldnames]
        writer = csv.DictWriter(output, fieldnames=fields + ["selected_primary_term"], delimiter="\t")
        writer.writeheader()
        for row in reader:
            display = selected_rows.pop(tuple(row[key] for key in keys), None)
            row["selected_primary_term"] = display is not None
            if display is not None:
                row.update({field: display[field] for field in selected_fields if field not in row})
            writer.writerow(row)
        for row in selected_rows.values():
            writer.writerow(dict(row, selected_primary_term=True))


def export_sensitivity_tables(root=REPO, output=None):
    root = Path(root).resolve()
    output = Path(output).resolve() if output is not None else root / "source_data/supplementary_tables/multilayer_sensitivity"
    required = [name for name, _, _ in SPECIFICATIONS] + [name for name, _ in PRIMARY_TERM_SELECTIONS.values()]
    missing = [name for name in required if not (root / name).is_file()]
    if missing:
        raise FileNotFoundError("Missing supplementary sources: " + ", ".join(missing))
    output.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(prefix=".sensitivity-tables-", dir=output.parent) as temporary:
        candidate = Path(temporary) / "tables"
        candidate.mkdir()
        index_rows = []
        for source_name, destination_name, figure in SPECIFICATIONS:
            source = root / source_name
            destination = candidate / destination_name
            destination.parent.mkdir(parents=True, exist_ok=True)
            if source_name in PRIMARY_TERM_SELECTIONS:
                selection, keys = PRIMARY_TERM_SELECTIONS[source_name]
                export_go_tests(source, root / selection, destination, keys)
            else:
                shutil.copyfile(source, destination)
            rows, columns = "", ""
            if ".tsv" in destination.suffixes:
                opener = gzip.open if destination.suffix == ".gz" else open
                with opener(destination, "rt") as handle:
                    reader = csv.reader(handle, delimiter="\t")
                    columns = len(next(reader))
                    rows = sum(1 for _ in reader)
            index_rows.append({"Figure": figure, "File": destination_name, "Rows": rows,
                               "Columns": columns})
        index_path = candidate / "table_index.tsv"
        with index_path.open("w") as handle:
            writer = csv.DictWriter(handle, fieldnames=["Figure", "File", "Rows", "Columns"], delimiter="\t")
            writer.writeheader()
            writer.writerows(index_rows)
        publish_directory(candidate, output)
    print(f"Exported {len(SPECIFICATIONS)} supplementary sensitivity tables.")
    return output


if __name__ == "__main__":
    export_sensitivity_tables()
