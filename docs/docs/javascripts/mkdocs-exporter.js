/**
 * An interface with the MkDocs Exporter plugin.
 *
 * To render mathjax when exporting page to PDFs.
 * https://github.com/adrienbrignon/mkdocs-exporter/blob/master/docs/assets/scripts/mkdocs-exporter.js
 */
window.MkDocsExporter = {

  /**
   * Render the page...
   */
  render: async () => {
    if (window.MathJax) {
      if (typeof window.MathJax.typesetPromise === 'function') {
        await window.MathJax.typesetPromise();
      }
    }
  }

};
