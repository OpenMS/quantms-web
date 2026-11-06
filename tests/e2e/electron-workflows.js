// quantms-web scenario for electron/scripts/e2e.js: runs the label-free DDA
// workflow (DecoyDatabase -> ProSE -> Percolator -> IDFilter -> ProteomicsLFQ)
// on two of the bundled example files, one per sample group, inside the
// Electron window, then checks that the Viewer draws their tables.
const { checkTable } = require("./viewer");

const RUN_TIMEOUT = 60 * 60 * 1000;
const FILES = [
  ["B000248_Ap_6883_EXT-765_DIA_Yeast_UPS2_ratio01_DDA.mzML"],
  ["B000264_Ap_6883_EXT-765_DIA_Yeast_UPS2_ratio10_DDA.mzML"],
];

module.exports = async (t) => {
  // Describe experiment fills a fresh workspace with the example data.
  await t.open("workflow_configure");
  // A fresh workspace selects every example file; keep two, one per group.
  // The demo files get their groups from their names.
  const files = t.page
    .locator('[data-testid="stMultiSelect"]')
    .filter({ has: t.page.locator('[data-testid="stWidgetLabel"]', { hasText: "MS data" }) })
    .first();
  await files.locator('[title="Clear all"], [aria-label="Clear all"]').first().click();
  await t.settle();
  for (const [file] of FILES) await t.choose("MS data", file);
  await t.shot("configured");

  await t.open("workflow_run");
  // The run page reruns every second while the workflow runs, so it
  // settles only when the run has ended.
  await t.clickButton("Start Workflow", { timeout: RUN_TIMEOUT });
  await t.waitForText("Workflow completed successfully.", {
    timeout: RUN_TIMEOUT,
    failOn: ["Errors occurred, check log file."],
  });
  await t.shot("workflow-finished");

  // The Viewer opens on the run that just finished, in its default layout.
  await t.open("viewer", { timeout: 10 * 60 * 1000 });
  await t.waitForComponent(".tabulator-row", { timeout: 5 * 60 * 1000 });
  await t.shot("viewer-default");

  await checkTable(t, "psms_search");
  await checkTable(t, "psms");
  await checkTable(t, "protein_abundance");
};
