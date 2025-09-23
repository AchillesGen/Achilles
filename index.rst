.. Achilles documentation master file, created by
   sphinx-quickstart on Mon Sep 22 23:17:23 2025.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to Achilles!
====================

.. raw:: html

   <div style="display: flex; flex-direction: column; align-items: center; margin-top: 2em;">

   <a href="latest/index.html" style="
       text-decoration: none;
       padding: 1em 2em;
       background-color: #4CAF50;
       color: white;
       border-radius: 8px;
       font-size: 1.2em;
       margin-bottom: 1em;
       display: inline-block;">
       📘 Go to Documentation
   </a>

   <label for="version-select" style="font-size: 1.1em; margin-bottom: 0.5em;">Choose a version:</label>
   <select id="version-select" onchange="if (this.value) window.location.href=this.value" style="
       font-size: 1em;
       padding: 0.5em;
       border-radius: 4px;">
   </select>

   <script>
     fetch('versions.json')
       .then(response => response.json())
       .then(versions => {
         const select = document.getElementById('version-select');
         versions.forEach(v => {
           const opt = document.createElement('option');
           opt.value = v + '/index.html';
           opt.textContent = v;
           if (v === 'latest') {
             opt.textContent = 'Latest Release';
             opt.selected = true;
           }
           select.appendChild(opt);
         });
       });
   </script>

   </div>
