import e, { Suspense as h, useEffect as k, useState as d } from "react";
import { o as g, L as _, u as y, r as N, i as I, j as M, n as R, I as j, M as B, t as p, av as L, aw as T, ax as q, ay as w, B as S, q as V } from "./index-8ee3baa4.js";
import "react-dom";
const W = g(function({ connectionType: a, model: o, session: t }) {
  const i = a.configEditorComponent || _;
  return e.createElement(
    h,
    { fallback: e.createElement(y.LoadingEllipses, null) },
    e.createElement(i, { model: { target: o }, session: t })
  );
});
var v = {}, z = I;
Object.defineProperty(v, "__esModule", {
  value: !0
});
var b = v.default = void 0, A = z(N()), D = M, O = (0, A.default)(/* @__PURE__ */ (0, D.jsx)("path", {
  d: "M19 19H5V5h7V3H5c-1.11 0-2 .9-2 2v14c0 1.1.89 2 2 2h14c1.1 0 2-.9 2-2v-7h-2v7zM14 3v2h3.59l-9.83 9.83 1.41 1.41L19 6.41V10h2V3h-7z"
}), "OpenInNew");
b = v.default = O;
const F = g(function({ connectionTypeChoices: o, connectionType: t, setConnectionType: i }) {
  return k(() => {
    t || i(o[0]);
  }), e.createElement("form", { autoComplete: "off" }, t ? e.createElement(R, { value: t.name, label: "connectionType", helperText: t.description ? e.createElement(
    e.Fragment,
    null,
    t.description,
    t.url ? e.createElement(
      j,
      { href: t.url, rel: "noopener noreferrer", target: "_blank" },
      e.createElement(b, null)
    ) : null
  ) : null, select: !0, fullWidth: !0, onChange: (n) => i(o.find((u) => u.name === n.target.value)), variant: "outlined" }, o.map((n) => e.createElement(B, { key: n.name, value: n.name }, n.displayName || n.name))) : null);
}), H = V()((a) => ({
  root: {
    marginTop: a.spacing(1)
  },
  stepper: {
    backgroundColor: a.palette.background.default
  },
  button: {
    marginTop: a.spacing(1),
    marginRight: a.spacing(1)
  },
  actionsContainer: {
    marginBottom: a.spacing(2)
  }
})), f = ["Select a Connection Type", "Configure Connection"], J = g(function({ model: o }) {
  const [t, i] = d(), [n, u] = d(), [r, C] = d(0), { classes: c } = H(), s = p.getSession(o), { pluginManager: x } = p.getEnv(s);
  return e.createElement(
    "div",
    { className: c.root },
    e.createElement(L, { className: c.stepper, activeStep: r, orientation: "vertical" }, f.map((E) => e.createElement(
      T,
      { key: E },
      e.createElement(q, null, E),
      e.createElement(
        w,
        null,
        r === 0 ? e.createElement(F, { connectionTypeChoices: x.getConnectionElements(), connectionType: t, setConnectionType: (l) => {
          if (i(l), !l)
            return;
          const m = `${l.name}-${Date.now()}`;
          u(l.configSchema.create({ connectionId: m }, p.getEnv(o)));
        } }) : t && n ? e.createElement(W, { connectionType: t, model: n, session: s }) : null,
        e.createElement(
          "div",
          { className: c.actionsContainer },
          e.createElement(S, { disabled: r === 0, onClick: () => C(r - 1), className: c.button }, "Back"),
          e.createElement(S, { disabled: !(r === 0 && t || r === 1 && n), variant: "contained", color: "primary", onClick: () => {
            var l;
            if (r === f.length - 1) {
              const m = s.addConnectionConf(n);
              (l = s.makeConnection) === null || l === void 0 || l.call(s, m), s.hideWidget(o);
            } else
              C(r + 1);
          }, className: c.button, "data-testid": "addConnectionNext" }, r === f.length - 1 ? "Connect" : "Next")
        )
      )
    )))
  );
});
export {
  J as default
};
