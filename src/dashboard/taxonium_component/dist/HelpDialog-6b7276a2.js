import e from "react";
import { u as l, b as n, J as r, c as a, B as o } from "./index-8ee3baa4.js";
import "react-dom";
function u({ handleClose: t }) {
  return e.createElement(
    l.Dialog,
    { open: !0, maxWidth: "xl", onClose: t, title: "Using the search box" },
    e.createElement(
      n,
      null,
      e.createElement("h3", null, "Searching"),
      e.createElement(
        "ul",
        null,
        e.createElement("li", null, "Jump to a feature or reference sequence by typing its name in the location box and pressing Enter."),
        e.createElement(
          "li",
          null,
          "Jump to a specific region by typing the region into the location box as: ",
          e.createElement("code", null, "ref:start..end"),
          " or ",
          e.createElement("code", null, "ref:start-end"),
          ". Commas are allowed in the start and end coordinates. A space-separated list of locstrings can be used to open up multiple chromosomes at a time"
        )
      ),
      e.createElement("h3", null, "Example Searches"),
      e.createElement(
        "ul",
        null,
        e.createElement(
          "li",
          null,
          e.createElement("code", null, "BRCA"),
          " - searches for the feature named BRCA"
        ),
        e.createElement(
          "li",
          null,
          e.createElement("code", null, "chr4"),
          " - jumps to chromosome 4"
        ),
        e.createElement(
          "li",
          null,
          e.createElement("code", null, "chr4:79,500,000..80,000,000"),
          " - jumps the region on chromosome 4 between 79.5Mb and 80Mb."
        ),
        e.createElement(
          "li",
          null,
          e.createElement("code", null, "chr1:1-100 chr2:1-100"),
          " - create a split view of chr1:1-100 and chr2:1-100"
        ),
        e.createElement(
          "li",
          null,
          e.createElement("code", null, "chr1 chr2 chr3"),
          " - open up multiple chromosomes at once"
        ),
        e.createElement(
          "li",
          null,
          e.createElement("code", null, "chr1:1-100[rev] chr2:1-100"),
          " - open up the first region in the horizontally flipped orientation"
        ),
        e.createElement(
          "li",
          null,
          e.createElement("code", null, "chr1 100 200"),
          " - use whitespace separated refname, start, end"
        )
      )
    ),
    e.createElement(r, null),
    e.createElement(
      a,
      null,
      e.createElement(o, { onClick: () => t(), color: "primary" }, "Close")
    )
  );
}
export {
  u as default
};
