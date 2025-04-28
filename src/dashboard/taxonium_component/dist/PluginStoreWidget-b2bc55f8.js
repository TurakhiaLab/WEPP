import e, { useState as v, useEffect as K } from "react";
import { r as x, i as S, j as _, q as $, o as b, t as h, aB as F, aC as O, a4 as Q, I, U as z, T as g, a9 as X, D as ee, a as te, b as B, c as q, B as P, a3 as ae, ab as ne, ad as le, Z as re, aD as se, p as oe, H as ie, u as W, aE as j, n as y, an as D, aF as ce, a8 as ue, a5 as me, al as T, am as A } from "./index-8ee3baa4.js";
import { d as de } from "./Clear-70a3606f.js";
import "react-dom";
var w = {}, ge = S;
Object.defineProperty(w, "__esModule", {
  value: !0
});
var J = w.default = void 0, fe = ge(x()), ve = _, pe = (0, fe.default)(/* @__PURE__ */ (0, ve.jsx)("path", {
  d: "M11 7h2v2h-2zm0 4h2v6h-2zm1-9C6.48 2 2 6.48 2 12s4.48 10 10 10 10-4.48 10-10S17.52 2 12 2zm0 18c-4.41 0-8-3.59-8-8s3.59-8 8-8 8 3.59 8 8-3.59 8-8 8z"
}), "InfoOutlined");
J = w.default = pe;
var M = {}, Ee = S;
Object.defineProperty(M, "__esModule", {
  value: !0
});
var H = M.default = void 0, he = Ee(x()), Pe = _, Ce = (0, he.default)(/* @__PURE__ */ (0, Pe.jsx)("path", {
  d: "M18 8h-1V6c0-2.76-2.24-5-5-5S7 3.24 7 6v2H6c-1.1 0-2 .9-2 2v10c0 1.1.9 2 2 2h12c1.1 0 2-.9 2-2V10c0-1.1-.9-2-2-2zm-6 9c-1.1 0-2-.9-2-2s.9-2 2-2 2 .9 2 2-.9 2-2 2zm3.1-9H8.9V6c0-1.71 1.39-3.1 3.1-3.1 1.71 0 3.1 1.39 3.1 3.1v2z"
}), "Lock");
H = M.default = Ce;
const V = $()(() => ({
  closeDialog: {
    position: "absolute",
    right: 0,
    top: 0
  },
  dialogContainer: {
    margin: 15
  },
  lockedPluginTooltip: {
    marginRight: "0.5rem"
  }
}));
function ye() {
  const { classes: t } = V();
  return e.createElement(
    X,
    { className: t.lockedPluginTooltip, title: "This plugin was installed by an administrator, you cannot remove it." },
    e.createElement(H, null)
  );
}
function be({ onClose: t, plugin: a }) {
  const { classes: r } = V();
  return e.createElement(
    ee,
    { open: !0, onClose: () => t() },
    e.createElement(
      te,
      null,
      e.createElement(
        I,
        { className: r.closeDialog, "aria-label": "close-dialog", onClick: () => t() },
        e.createElement(z, null)
      )
    ),
    e.createElement(
      B,
      null,
      e.createElement(
        g,
        null,
        "Please confirm that you want to remove ",
        a,
        ". Note: if any resources in this session still use this plugin, it may cause your session to crash"
      ),
      e.createElement(
        q,
        null,
        e.createElement(P, { variant: "contained", color: "primary", onClick: () => {
          window.setTimeout(() => {
            t(a);
          }, 500);
        } }, "Confirm"),
        e.createElement(P, { variant: "contained", color: "secondary", onClick: () => {
          t();
        } }, "Cancel")
      )
    )
  );
}
function xe({ plugin: t, model: a, pluginManager: r }) {
  const [c, o] = v(), n = h.getSession(a), { sessionPlugins: i } = n, m = i == null ? void 0 : i.some((u) => r.pluginMetadata[t.name].url === u.url), d = F(a, 3), { jbrowse: f, adminMode: l } = d;
  return e.createElement(
    e.Fragment,
    null,
    c ? e.createElement(be, { plugin: c, onClose: (u) => {
      if (u) {
        const s = r.pluginMetadata[t.name];
        l ? f.removePlugin(s) : O.isSessionWithSessionPlugins(n) && n.removeSessionPlugin(s);
      }
      o(void 0);
    } }) : null,
    e.createElement(
      Q,
      { key: t.name },
      l || m ? e.createElement(
        I,
        { "aria-label": "removePlugin", "data-testid": `removePlugin-${t.name}`, onClick: () => o(t.name) },
        e.createElement(z, null)
      ) : e.createElement(ye, null),
      e.createElement(g, null, t.name)
    )
  );
}
const Se = b(xe);
function _e({ pluginManager: t, model: a }) {
  const { plugins: r } = t, c = new Set(r.filter((n) => {
    var i;
    return (i = t.pluginMetadata[n.name]) === null || i === void 0 ? void 0 : i.isCore;
  }).map((n) => n.name)), o = r.filter((n) => !c.has(n.name));
  return e.createElement(ae, null, o.length > 0 ? o.filter((n) => n.name.toLowerCase().includes(a.filterText.toLowerCase())).map((n) => e.createElement(Se, { key: n.name, plugin: n, model: a, pluginManager: t })) : e.createElement(g, null, "No plugins currently installed"));
}
const $e = b(_e);
var R = {}, Ie = S;
Object.defineProperty(R, "__esModule", {
  value: !0
});
var Y = R.default = void 0, ke = Ie(x()), je = _, De = (0, ke.default)(/* @__PURE__ */ (0, je.jsx)("path", {
  d: "M12 12c2.21 0 4-1.79 4-4s-1.79-4-4-4-4 1.79-4 4 1.79 4 4 4zm0 2c-2.67 0-8 1.34-8 4v2h16v-2c0-2.66-5.33-4-8-4z"
}), "Person");
Y = R.default = De;
var U = {}, we = S;
Object.defineProperty(U, "__esModule", {
  value: !0
});
var Z = U.default = void 0, Me = we(x()), Re = _, Ue = (0, Me.default)(/* @__PURE__ */ (0, Re.jsx)("path", {
  d: "M9 16.17 4.83 12l-1.42 1.41L9 19 21 7l-1.41-1.41z"
}), "Check");
Z = U.default = Ue;
const Le = $()({
  card: {
    margin: "1em"
  },
  icon: {
    marginLeft: "0.5em",
    marginRight: "0.5em"
  },
  bold: {
    fontWeight: 600
  },
  dataField: {
    display: "flex",
    alignItems: "center",
    margin: "0.4em 0em"
  }
}), Ne = b(function({ plugin: a, model: r, adminMode: c }) {
  const { classes: o } = Le(), n = h.getSession(r), { pluginManager: i } = h.getEnv(r), { runtimePluginDefinitions: m } = i, d = m.some((C) => C.url === a.url), [f, l] = v(!1), u = d || f, s = F(r, 3), { jbrowse: p } = s;
  return e.createElement(
    ne,
    { variant: "outlined", key: a.name, className: o.card },
    e.createElement(
      le,
      null,
      e.createElement(
        "div",
        { className: o.dataField },
        e.createElement(
          g,
          { variant: "h5" },
          e.createElement(re, { href: `${a.location}#readme`, target: "_blank", rel: "noopener" }, a.name)
        )
      ),
      e.createElement(
        "div",
        { className: o.dataField },
        e.createElement(Y, { style: { marginRight: "0.5em" } }),
        e.createElement(g, null, a.authors.join(", "))
      ),
      e.createElement(g, { className: o.bold }, "Description:"),
      e.createElement(g, null, a.description)
    ),
    e.createElement(
      se,
      null,
      e.createElement(P, { variant: "contained", disabled: u, startIcon: d ? e.createElement(Z, null) : e.createElement(oe, null), onClick: () => {
        c ? p.addPlugin({ name: a.name, url: a.url }) : O.isSessionWithSessionPlugins(n) && n.addSessionPlugin(a), l(!0);
      } }, d ? "Installed" : "Install")
    )
  );
}), Te = $()((t) => ({
  dialogContent: {
    display: "flex",
    flexDirection: "column"
  },
  expand: {
    transform: "rotate(0deg)",
    marginLeft: "auto",
    transition: t.transitions.create("transform", {
      duration: t.transitions.duration.shortest
    })
  },
  expandOpen: {
    transform: "rotate(180deg)"
  }
}));
function Ae({ open: t, onClose: a, model: r }) {
  const { classes: c, cx: o } = Te(), [n, i] = v(""), [m, d] = v(""), [f, l] = v(""), [u, s] = v(""), [p, C] = v(!1), { jbrowse: k } = ie(r), G = !!(n && m || f || u);
  function L() {
    n && m && k.addPlugin({ name: n, umdUrl: m }), f && k.addPlugin({ esmUrl: f }), u && k.addPlugin({ cjsUrl: u });
  }
  function N() {
    i(""), d(""), l(""), s(""), a();
  }
  return e.createElement(
    W.Dialog,
    { open: t, onClose: N, title: "Add custom plugin" },
    e.createElement(
      "form",
      { onSubmit: L },
      e.createElement(
        B,
        { className: c.dialogContent },
        e.createElement(j, null, "Enter the name of the plugin and its URL. The name should match what is defined in the plugin's build."),
        e.createElement(y, { id: "umd-name-input", name: "umdName", label: "Plugin name", variant: "outlined", value: n, onChange: (E) => i(E.target.value) }),
        e.createElement(y, { id: "umd-url-input", name: "umdUrl", label: "Plugin URL", variant: "outlined", value: m, onChange: (E) => d(E.target.value) }),
        e.createElement(
          j,
          { onClick: () => C(!p) },
          e.createElement(
            I,
            { className: o(c.expand, {
              [c.expandOpen]: p
            }), "aria-expanded": p, "aria-label": "show more" },
            e.createElement(D, null)
          ),
          "Advanced options"
        ),
        e.createElement(
          ce,
          { in: p },
          e.createElement(
            "div",
            { className: c.dialogContent },
            e.createElement(j, null, "The above fields assume that the plugin is built in UMD format. If your plugin is in another format, or you have additional builds you want to add (such as a CJS build for using NodeJS APIs in desktop), you can enter the URLs for those builds below."),
            e.createElement(y, { id: "esm-url-input", name: "esmUrl", label: "ESM build URL", variant: "outlined", value: f, onChange: (E) => l(E.target.value) }),
            e.createElement(y, { id: "cjs-url-input", name: "cjsUrl", label: "CJS build URL", variant: "outlined", value: u, onChange: (E) => s(E.target.value) })
          )
        )
      ),
      e.createElement(
        q,
        null,
        e.createElement(P, { variant: "contained", onClick: N }, "Cancel"),
        e.createElement(P, { variant: "contained", color: "primary", onClick: L, disabled: !G }, "Submit")
      )
    )
  );
}
const Fe = b(Ae), Oe = $()((t) => ({
  root: {
    margin: t.spacing(1)
  },
  expandIcon: {
    color: t.palette.tertiary.contrastText
  },
  adminBadge: {
    margin: "0.5em",
    borderRadius: 3,
    backgroundColor: t.palette.quaternary.main,
    padding: "1em",
    display: "flex",
    alignContent: "center"
  },
  customPluginButton: {
    margin: "0.5em",
    display: "flex",
    justifyContent: "center"
  }
}));
function ze({ model: t }) {
  const { classes: a } = Oe(), [r, c] = v(), [o, n] = v(), [i, m] = v(!1), { adminMode: d } = h.getSession(t), { pluginManager: f } = ue(t);
  return K(() => {
    const l = new AbortController(), { signal: u } = l;
    return (async () => {
      try {
        const s = await fetch("https://jbrowse.org/plugin-store/plugins.json", { signal: u });
        if (!s.ok) {
          const C = await s.text();
          throw new Error(`Failed to fetch plugin data: ${s.status} ${s.statusText} ${C}`);
        }
        const p = await s.json();
        u.aborted || c(p.plugins);
      } catch (s) {
        console.error(s), n(s);
      }
    })(), () => {
      l.abort();
    };
  }, []), e.createElement(
    "div",
    { className: a.root },
    d && e.createElement(
      e.Fragment,
      null,
      !h.isElectron && e.createElement(
        "div",
        { className: a.adminBadge },
        e.createElement(J, { style: { marginRight: "0.3em" } }),
        e.createElement(
          g,
          null,
          "You are using the ",
          e.createElement("code", null, "admin-server"),
          ". Any changes you make will be saved to your configuration file. You also have the ability to add custom plugins that are not in the store."
        )
      ),
      e.createElement(
        "div",
        { className: a.customPluginButton },
        e.createElement(P, { variant: "contained", onClick: () => m(!0) }, "Add custom plugin")
      ),
      e.createElement(Fe, { open: i, onClose: () => m(!1), model: t })
    ),
    e.createElement(y, { label: "Filter plugins", value: t.filterText, onChange: (l) => t.setFilterText(l.target.value), fullWidth: !0, InputProps: {
      endAdornment: e.createElement(
        me,
        { position: "end" },
        e.createElement(
          I,
          { onClick: () => t.clearFilterText() },
          e.createElement(de, null)
        )
      )
    } }),
    e.createElement(
      T,
      { defaultExpanded: !0 },
      e.createElement(
        A,
        { expandIcon: e.createElement(D, { className: a.expandIcon }) },
        e.createElement(g, { variant: "h5" }, "Installed plugins")
      ),
      e.createElement(
        "div",
        { style: { margin: "1em" } },
        e.createElement($e, { pluginManager: f, model: t })
      )
    ),
    e.createElement(
      T,
      { defaultExpanded: !0 },
      e.createElement(
        A,
        { expandIcon: e.createElement(D, { className: a.expandIcon }) },
        e.createElement(g, { variant: "h5" }, "Available plugins")
      ),
      o ? e.createElement(g, { color: "error" }, `${o}`) : r ? r.filter((l) => !(h.isElectron && l.cjsUrl) && l.name.toLowerCase().includes(t.filterText.toLowerCase())).map((l) => e.createElement(Ne, { key: l.name, plugin: l, model: t, adminMode: !!d })) : e.createElement(W.LoadingEllipses, null)
    )
  );
}
const He = b(ze);
export {
  He as default
};
