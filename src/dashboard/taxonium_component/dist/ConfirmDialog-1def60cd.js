import e, { useState as u } from "react";
import { _ as m, b as s, T as r, l as f, n as p, c as g, B as o } from "./index-8ee3baa4.js";
import "react-dom";
const d = ({ tracks: l, onClose: a }) => {
  const [n, i] = u(`MultiWiggle ${Date.now()}`), c = l.every((t) => t.type === "QuantitativeTrack");
  return e.createElement(
    m,
    { open: !0, onClose: () => a(!1), title: "Confirm multi-wiggle track create" },
    e.createElement(
      s,
      null,
      e.createElement(
        r,
        null,
        c ? null : "Not every track looks like a QuantitativeTrack. This could have unexpected behavior, confirm if it looks ok.",
        "Listing:"
      ),
      e.createElement("ul", null, l.map((t) => e.createElement(
        "li",
        { key: t.trackId },
        f.readConfObject(t, "name"),
        " - ",
        t.type
      ))),
      e.createElement(p, { value: n, onChange: (t) => i(t.target.value), helperText: "Track name" }),
      e.createElement(r, null, "Confirm creation of track?")
    ),
    e.createElement(
      g,
      null,
      e.createElement(o, { onClick: () => a(!1), color: "primary" }, "Cancel"),
      e.createElement(o, { onClick: () => a(!0, { name: n }), color: "primary", variant: "contained", autoFocus: !0 }, "Submit")
    )
  );
};
export {
  d as default
};
