import e, { useState as d } from "react";
import { o as m, u as g, b as p, T as r, n as u, c as f, B as l } from "./index-8ee3baa4.js";
import "react-dom";
function h(n) {
  const { model: c, handleClose: a } = n, [t, s] = d(""), o = t.match(/^[A-Za-z][A-Za-z0-9]$/);
  return e.createElement(
    g.Dialog,
    { open: !0, onClose: a, title: "Color by tag" },
    e.createElement(
      p,
      { style: { overflowX: "hidden" } },
      e.createElement(r, null, "Enter tag to color by: "),
      e.createElement(r, { color: "textSecondary" }, "Examples: XS or TS for RNA-seq inferred read strand, ts (lower-case) for minimap2 read strand, HP for haplotype, RG for read group, etc."),
      e.createElement(u, { value: t, onChange: (i) => s(i.target.value), placeholder: "Enter tag name", inputProps: {
        maxLength: 2,
        "data-testid": "color-tag-name-input"
      }, error: t.length === 2 && !o, helperText: t.length === 2 && !o ? "Not a valid tag" : "", autoComplete: "off", "data-testid": "color-tag-name" }),
      e.createElement(
        f,
        null,
        e.createElement(l, { variant: "contained", color: "primary", onClick: () => {
          c.setColorScheme({
            type: "tag",
            tag: t
          }), a();
        }, disabled: !o }, "Submit"),
        e.createElement(l, { variant: "contained", color: "secondary", onClick: a }, "Cancel")
      )
    )
  );
}
const T = m(h);
export {
  T as default
};
