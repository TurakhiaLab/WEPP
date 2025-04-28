import t, { useState as m } from "react";
import { o as g, u, b as p, T as r, n as d, c as y, B as n } from "./index-8ee3baa4.js";
import "react-dom";
function h(l) {
  const { model: c, handleClose: a } = l, [e, s] = m(""), o = e.match(/^[A-Za-z][A-Za-z0-9]$/);
  return t.createElement(
    u.Dialog,
    { open: !0, onClose: a, title: "Sort by tag" },
    t.createElement(
      p,
      null,
      t.createElement(r, null, "Set the tag to sort by"),
      t.createElement(r, { color: "textSecondary" }, "Examples: HP for haplotype, RG for read group, etc."),
      t.createElement(d, { value: e, onChange: (i) => s(i.target.value), placeholder: "Enter tag name", inputProps: {
        maxLength: 2,
        "data-testid": "sort-tag-name-input"
      }, error: e.length === 2 && !o, helperText: e.length === 2 && !o ? "Not a valid tag" : "", autoComplete: "off", "data-testid": "sort-tag-name" }),
      t.createElement(
        y,
        null,
        t.createElement(n, { variant: "contained", color: "primary", type: "submit", autoFocus: !0, onClick: () => {
          c.setSortedBy("tag", e), a();
        } }, "Submit"),
        t.createElement(n, { variant: "contained", color: "secondary", onClick: () => a() }, "Cancel")
      )
    )
  );
}
const S = g(h);
export {
  S as default
};
