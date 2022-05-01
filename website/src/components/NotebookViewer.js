import React from 'react';
import NbViewer from 'react-nbviewer';
import axios from "axios";
import ReactMarkdown from 'react-markdown';
import remarkMath from 'remark-math';
import rehypeKatex from 'rehype-katex';
import SyntaxHighlighter from 'react-syntax-highlighter';

import 'katex/dist/katex.min.css';
import 'react-nbviewer/dist/index.css';

// Cleans up json to remove newlines from math envs
const clean_up_json = (object) => {

  let metadata = object.metadata;

  const language =
    metadata.kernelspec.language == "Wolfram Language"
      ? "mathematica"
      : metadata.kernelspec.language;

  metadata.kernelspec.language = language;

  return {
    cells: object.cells.map((el) => {
      if (el.cell_type == "markdown") {
        const source = el.source.join("");
        const clean = source
          .trim()
          .replace(
            /(\${1,2})((?:\\.|[\s\S])*?)\1/g,
            (m, tag, src) => tag + src.replace(/\r?\n/g, "") + tag
          );

        return {
          cell_type: el.cell_type,
          id: el.id,
          metadata: el.metadata,
          source: clean.split(/(\r\n|\n|\r)/)
        };
      } else {
        return el;
      }
    }),
    metadata: metadata,
    nbformat: object.nbformat,
    nbformat_minor: object.nbformat_minor
  };
};

// Decode unicode characters
const b64DecodeUnicode = (str) => {
  return decodeURIComponent(
    atob(str)
      .split("")
      .map(function (c) {
        return "%" + ("00" + c.charCodeAt(0).toString(16)).slice(-2);
      })
      .join("")
  );
};

// Read clean JSON data from git .ipynb url
const getNotebook = async (notebookURL) => {
  const regex = new RegExp(
    /https:\/\/github.com\/(.+?)\/(.+?)\/blob\/(.+?)\/(.+\.ipynb)/
  );
  const parts = notebookURL.match(regex);
  const [, organization, repo, branch, notebook] = parts;
  const requestURL = `https://api.github.com/repos/${organization}/${repo}/contents/${notebook}?ref=${branch}`;
  const data = await (await axios.get(requestURL)).data;
  const content = clean_up_json(JSON.parse(b64DecodeUnicode(data.content)));
  //const content = JSON.parse(b64DecodeUnicode(data.content)); // disable TeX rendering
  console.log(content);
  return content;
};

// Markdown Components
//
// TeX Support
const MathMarkdown = (props) => <ReactMarkdown 
	remarkPlugins={[remarkMath]} 
	rehypePlugins={[rehypeKatex]} 
	children={props.source} 
	/>
	
// source -> children breaking change since 5.0.3
const Markdown = (props) => <ReactMarkdown children={props.source} />

// Notebook Viewer
export default function NotebookViewer({ notebookURL }) {
  const [notebook, setNotebook] = React.useState(null);

  React.useEffect(() => {
    getNotebook(notebookURL).then(setNotebook);
  }, [notebookURL]);

  return (
    <div className="NotebookViewer">
      {notebook === null ? (
        "Loading notebook..."
      ) : (
        <NbViewer
          source={notebook}
          code={SyntaxHighlighter}
          markdown={MathMarkdown}
          //markdown={Markdown} // disable TeX rendering
        />
      )}
    </div>
  );
}

