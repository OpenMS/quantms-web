// quantms-web scenario for electron/scripts/e2e.js: runs the label-free DDA
// workflow (DecoyDatabase -> Comet -> Percolator -> IDFilter -> ProteomicsLFQ)
// on two of the bundled example files, one per sample group, inside the
// Electron window, then checks that the result pages draw their tables.
const RUN_TIMEOUT = 60 * 60 * 1000;
const FILES = [
  ["B000248_Ap_6883_EXT-765_DIA_Yeast_UPS2_ratio01_DDA.mzML", "ratio01"],
  ["B000264_Ap_6883_EXT-765_DIA_Yeast_UPS2_ratio10_DDA.mzML", "ratio10"],
];

async function fillText(t, label, value) {
  const box = t.page.getByRole("textbox", { name: label, exact: true }).first();
  await box.fill(value);
  await box.press("Enter");
  await t.settle();
}

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
  for (const [file] of FILES) await t.choose("mzML-files", file);
  await t.clickTab("Group Selection");
  for (const [file, group] of FILES) await fillText(t, `Group for ${file}`, group);
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
