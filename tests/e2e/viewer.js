// Helpers for the quantms-web scenarios: the Viewer shows the components of
// the saved layout, so a scenario imports a layout through Edit layout to see
// the component it checks.
const fs = require("fs");
const path = require("path");

// Opens the Viewer showing `rows` (rows of component names) for the run
// picked above it, which after a run is the run that just finished.
async function showComponents(t, rows, name) {
  await t.open("viewer", { timeout: 10 * 60 * 1000 });
  await t.clickButton("Edit layout");
  const file = path.join(t.out, `${name}-layout.json`);
  fs.writeFileSync(file, JSON.stringify([rows]));
  const editor = await t.waitForComponent("#file");
  await editor.locator("#file").setInputFiles(file);
  await editor.getByRole("button", { name: "Save", exact: true }).click();
  // Save reruns the page, which closes the dialog and draws the new layout.
  await t.page.locator('[role="dialog"]').waitFor({ state: "detached", timeout: 120000 });
  await t.settle({ timeout: 10 * 60 * 1000 });
}

// Shows `component` alone and waits for its table to draw rows.
async function checkTable(t, component) {
  await showComponents(t, [[component]], component);
  // The first draw after a run builds the components' caches, which took
  // over two minutes on a macOS runner.
  const frame = await t.waitForComponent(".tabulator-row", { timeout: 5 * 60 * 1000 });
  const rows = await frame.locator(".tabulator-row").count();
  t.log(`${component}: ${rows} table rows`);
  await t.shot(`${component}-results`);
}

module.exports = { showComponents, checkTable };
