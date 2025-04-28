import e from "react";
import { u as s, aU as l, b as i, T as r, Z as t, q as u } from "./index-8ee3baa4.js";
import "react-dom";
const c = u()({
  content: {
    minWidth: 800
  }
});
function b({ open: n, onClose: o }) {
  const { classes: a } = c();
  return e.createElement(
    s.Dialog,
    { open: n, onClose: o, maxWidth: "xl", title: `JBrowse v${l}` },
    e.createElement(
      i,
      { className: a.content },
      e.createElement(r, null, "JBrowse is a GMOD project © 2019-2021, The Evolutionary Software Foundation"),
      e.createElement(r, null, "Here are some resources and documentation. Please report the version number above when asking questions. Thanks!"),
      e.createElement(
        "ul",
        null,
        e.createElement(
          "li",
          null,
          e.createElement(t, { href: "https://github.com/GMOD/jbrowse-components/discussions", target: "_blank", rel: "noopener noreferrer" }, "Question & answer forum")
        ),
        e.createElement(
          "li",
          null,
          e.createElement(t, { href: "https://github.com/GMOD/jbrowse-components/issues/new/choose", target: "_blank", rel: "noopener noreferrer" }, "Report a bug")
        ),
        e.createElement(
          "li",
          null,
          e.createElement(t, { href: "https://jbrowse.org/jb2/docs/user_guide", target: "_blank", rel: "noopener noreferrer" }, "User guide")
        ),
        e.createElement(
          "li",
          null,
          e.createElement(t, { href: "https://jbrowse.org/jb2/docs/", target: "_blank", rel: "noopener noreferrer" }, "Documentation")
        )
      )
    )
  );
}
export {
  b as default
};
