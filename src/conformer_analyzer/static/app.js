const input = document.querySelector("#log-files");
const list = document.querySelector("#file-list");
const dropZone = document.querySelector("#drop-zone");
const form = document.querySelector("#analysis-form");
let selectedFiles = [];

function syncInput() {
  const transfer = new DataTransfer();
  selectedFiles.forEach((file) => transfer.items.add(file));
  input.files = transfer.files;
}

function renderFiles() {
  list.replaceChildren();
  selectedFiles.forEach((file, index) => {
    const item = document.createElement("li");
    const name = document.createElement("span");
    name.textContent = file.name;
    item.appendChild(name);
    [["↑", -1, "Move up"], ["↓", 1, "Move down"], ["×", 0, "Remove"]].forEach(([symbol, offset, label]) => {
      const button = document.createElement("button");
      button.type = "button";
      button.textContent = symbol;
      button.title = label;
      button.disabled = (offset === -1 && index === 0) || (offset === 1 && index === selectedFiles.length - 1);
      button.addEventListener("click", () => {
        if (offset === 0) selectedFiles.splice(index, 1);
        else [selectedFiles[index], selectedFiles[index + offset]] = [selectedFiles[index + offset], selectedFiles[index]];
        syncInput();
        renderFiles();
      });
      item.appendChild(button);
    });
    list.appendChild(item);
  });
}

input.addEventListener("change", () => {
  selectedFiles = Array.from(input.files);
  renderFiles();
});

["dragenter", "dragover"].forEach((eventName) => dropZone.addEventListener(eventName, (event) => {
  event.preventDefault();
  dropZone.classList.add("dragging");
}));
["dragleave", "drop"].forEach((eventName) => dropZone.addEventListener(eventName, (event) => {
  event.preventDefault();
  dropZone.classList.remove("dragging");
}));
dropZone.addEventListener("drop", (event) => {
  selectedFiles = Array.from(event.dataTransfer.files).filter((file) => {
    const name = file.name.toLowerCase();
    return name.endsWith(".log") || name.endsWith(".xyz");
  });
  syncInput();
  renderFiles();
});

form.addEventListener("submit", () => {
  document.querySelector("#analyze-button").disabled = true;
  document.querySelector("#working").hidden = false;
});
