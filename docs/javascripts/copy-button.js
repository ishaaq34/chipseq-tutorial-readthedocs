// Add copy buttons to all code blocks
document.addEventListener('DOMContentLoaded', function() {
    // Find all code blocks
    const codeBlocks = document.querySelectorAll('pre code, pre');
    
    codeBlocks.forEach(function(block) {
        // Skip if already has a copy button
        if (block.parentElement.querySelector('.copy-button')) {
            return;
        }
        
        // Get the actual pre element
        const pre = block.tagName === 'PRE' ? block : block.parentElement;
        
        // Skip if not a pre element
        if (pre.tagName !== 'PRE') {
            return;
        }
        
        // Create copy button
        const button = document.createElement('button');
        button.className = 'copy-button';
        button.textContent = 'Copy';
        button.setAttribute('aria-label', 'Copy code to clipboard');
        
        // Add click handler
        button.addEventListener('click', async function() {
            // Get code text
            const code = pre.querySelector('code') || pre;
            const text = code.textContent;
            
            try {
                // Copy to clipboard
                await navigator.clipboard.writeText(text);
                
                // Show success state
                button.textContent = 'Copied!';
                button.classList.add('copied');
                
                // Reset after 2 seconds
                setTimeout(function() {
                    button.textContent = 'Copy';
                    button.classList.remove('copied');
                }, 2000);
            } catch (err) {
                // Fallback for older browsers
                const textarea = document.createElement('textarea');
                textarea.value = text;
                textarea.style.position = 'fixed';
                textarea.style.opacity = '0';
                document.body.appendChild(textarea);
                textarea.select();
                
                try {
                    document.execCommand('copy');
                    button.textContent = 'Copied!';
                    button.classList.add('copied');
                    
                    setTimeout(function() {
                        button.textContent = 'Copy';
                        button.classList.remove('copied');
                    }, 2000);
                } catch (e) {
                    button.textContent = 'Failed';
                    setTimeout(function() {
                        button.textContent = 'Copy';
                    }, 2000);
                }
                
                document.body.removeChild(textarea);
            }
        });
        
        // Add button to pre element
        pre.style.position = 'relative';
        pre.appendChild(button);
    });
});
