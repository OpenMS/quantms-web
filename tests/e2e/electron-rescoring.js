// quantms-web scenario for electron/scripts/e2e.js, run after
// electron-workflows.js (which leaves two example files selected): turns on
// "Add predicted-spectrum features", runs the label-free workflow again and
// checks that AlphaPeptDeep predicted the matches inside the packaged app.
// The first use downloads the AlphaPeptDeep models, so the runner needs
// network access.
const RUN_TIMEOUT = 60 * 60 * 1000;
const SWITCH = "Add predicted-spectrum features";

module.exports = async (t) => {
  await t.open("workflow_configure");
  await t.clickTab("3 · Rescoring");
  const body = await t.page.locator("body").innerText();
  if (body.includes("is not installed here")) {
    throw new Error("AlphaPeptDeep is not installed in this build");
  }
  const box = t.page
    .locator('[data-testid="stCheckbox"]')
    .filter({ hasText: SWITCH })
    .first();
  if (!(await box.locator("input").isChecked())) {
    await box.locator("label").first().click();
    await t.settle();
  }
  if (!(await box.locator("input").isChecked())) throw new Error(`"${SWITCH}" did not switch on`);
  await t.shot("rescoring-switch-on");

  await t.open("workflow_run");
  await t.clickButton("Start Workflow");
  await t.waitForText("Workflow completed successfully.", {
    timeout: RUN_TIMEOUT,
    failOn: ["Errors occurred, check log file."],
  });
  await t.shot("rescoring-run-finished");

  await t.open("results_rescoring", { timeout: 10 * 60 * 1000 });
  await t.waitForText("Rescored with predicted spectra and retention times.", { timeout: 5 * 60 * 1000 });
  // One summary row per file: "... | NCE 30 | 1,234 of 1,300 | 0.85 | 0.26 |"
  const text = await t.page.locator("body").innerText();
  const rows = [...text.matchAll(/NCE [\d.]+\s+([\d,]+) of ([\d,]+)\s+([\d.]+|None)\s+([\d.]+|None)/g)];
  if (rows.length !== 2) throw new Error(`expected 2 summary rows, found ${rows.length}`);
  for (const [, predicted, , targets, decoys] of rows) {
    t.log(`predicted ${predicted}, median spectral angle targets ${targets}, decoys ${decoys}`);
    if (Number(predicted.replace(/,/g, "")) === 0) throw new Error("no match was predicted");
    // Correct matches look like their prediction, decoys do not
    if (!(Number(targets) > Number(decoys))) throw new Error(`targets ${targets} not above decoys ${decoys}`);
  }
  await t.shot("rescoring-summary");
};
