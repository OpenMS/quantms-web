// quantms-web scenario for electron/scripts/e2e.js: runs the label-free DDA
// workflow (DecoyDatabase -> Comet -> Percolator -> IDFilter -> ProteomicsLFQ)
// on two of the bundled example files, one per sample group, inside the
// Electron window, then checks that the result pages draw their tables.
const RUN_TIMEOUT = 60 * 60 * 1000;
const FILES = [
  ["B000248_Ap_6883_EXT-765_DIA_Yeast_UPS2_ratio01_DDA.mzML"],
  ["B000264_Ap_6883_EXT-765_DIA_Yeast_UPS2_ratio10_DDA.mzML"],
];

async function checkTables(t, page) {
  await t.open(page);
  const frame = await t.waitForComponent(".tabulator-row", { timeout: 5 * 60 * 1000 });
  const rows = await frame.locator(".tabulator-row").count();
  t.log(`${page}: ${rows} table rows`);
  await t.shot(`${page}-results`);
}

module.exports = async (t) => {
  // The upload page fills the workspace with the example data.
  await t.open("workflow_fileupload");
  await t.open("workflow_configure");
  // A fresh workspace selects every example file; keep two, one per group.
  // The demo files get their groups from their names.
  const files = t.page
    .locator('[data-testid="stMultiSelect"]')
    .filter({ has: t.page.locator('[data-testid="stWidgetLabel"]', { hasText: "MS data files" }) })
    .first();
  await files.locator('[title="Clear all"], [aria-label="Clear all"]').first().click();
  await t.settle();
  for (const [file] of FILES) await t.choose("MS data files", file);
  await t.shot("configured");

  await t.open("workflow_run");
  await t.clickButton("Start Workflow");
  await t.waitForText("Workflow completed successfully.", {
    timeout: RUN_TIMEOUT,
    failOn: ["Errors occurred, check log file."],
  });
  await t.shot("workflow-finished");

  await checkTables(t, "results_database_search");
  await checkTables(t, "results_filtered");
  await checkTables(t, "results_abundance");
};
