// Helpers for the quantms-web scenarios: the Viewer shows the components of
// the current view, so a scenario makes a view of its own (+ in the view bar)
// holding the component it checks.

// Opens the Viewer showing `rows` (rows of component names, one component
// per row) for the run picked above it, which after a run is the run that
// just finished. `name` names the scenario's screenshots.
async function showComponents(t, rows, name) {
  await t.open("viewer", { timeout: 10 * 60 * 1000 });
  const bar = await t.waitForComponent("button.v.add");
  await bar.locator("button.v.add").click();
  // + opens Customize views on a copy of the current view: empty it, then
  // click each component in, every click adding a row
  const editor = await t.waitForComponent("#component-search");
  await editor.getByRole("button", { name: "Clear", exact: true }).click();
  for (const component of rows.flat()) {
    await editor.locator(`.chip[data-name="${component}"]`).click();
  }
  await t.shot(`${name}-view`);
  await editor.getByRole("button", { name: "Save", exact: true }).click();
  // Save reruns the page, which closes the dialog and draws the new view.
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
