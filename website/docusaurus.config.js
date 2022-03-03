const lightCodeTheme = require('prism-react-renderer/themes/github');
const darkCodeTheme = require('prism-react-renderer/themes/dracula');
const math = require('remark-math');
const katex = require('rehype-katex');

/** @type {import('@docusaurus/types').DocusaurusConfig} */
module.exports = {
  title: 'SpaRTaNS',
  tagline: 'Spatially Resolved Transport of Nonequilibrium Species',
  url: 'https://narang-lab.github.io',
  baseUrl: '/spartans/',
  onBrokenLinks: 'throw',
  onBrokenMarkdownLinks: 'warn',
  favicon: 'img/favicon.ico',
  organizationName: 'narang-lab', // Usually your GitHub org/user name.
  projectName: 'spartans', // Usually your repo name.
  trailingSlash: false, //gh-pages specific
  themeConfig: {
    navbar: {
      title: 'SpaRTaNS',
      logo: {
        alt: 'SpaRTaNS Logo',
        src: 'img/logo.svg',
      },
      items: [
        {
          type: 'doc',
          docId: 'intro',
          position: 'left',
          label: 'Docs',
        },
        {to: '/news', label: 'News', position: 'left'},
        {
          href: 'https://github.com/narang-lab/spartans',
          label: 'GitHub',
          position: 'right',
        },
      ],
    },
    footer: {
      style: 'dark',
      links: [
        {
          title: 'Docs',
          items: [
            {
              label: 'Documentation',
              to: '/docs/intro',
            },
          ],
        },
        {
          title: 'News',
          items: [
            {
              label: 'News',
              to: '/news',
            },
          ],
        },
        {
          title: 'More',
          items: [
            {
              label: 'Twitter',
              href: 'https://twitter.com/naranglab',
            },
            {
              label: 'GitHub',
              href: 'https://github.com/narang-lab/spartans',
            },
          ],
        },
      ],
      copyright: `Copyright © ${new Date().getFullYear()} NarangLab, Harvard University.`,
    },
    prism: {
      theme: lightCodeTheme,
      darkTheme: darkCodeTheme,
    },
  },
  presets: [
    [
      '@docusaurus/preset-classic',
      {
        docs: {
          sidebarPath: require.resolve('./sidebars.js'),
          // Please change this to your repo.
          editUrl: 'https://github.com/narang-lab/spartans/edit/main/website/',
          remarkPlugins: [math],
          rehypePlugins: [katex],
        },
        blog: {
          path: 'news',
          routeBasePath: 'news',
          showReadingTime: false,
          // Please change this to your repo.
          editUrl:
            'https://github.com/narang-lab/spartans/edit/main/website/news/',
        },
        theme: {
          customCss: require.resolve('./src/css/custom.css'),
        },
      },
    ],
  ],
};
