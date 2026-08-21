import { render } from "preact";
import { App } from "./app.tsx";
import "./index.css";

const root = document.getElementById("app");
if (!root) throw new Error("the #app mount point is missing");

render(<App />, root);
